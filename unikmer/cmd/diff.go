// Copyright © 2018-2021 Wei Shen <shenwei356@gmail.com>
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

package cmd

import (
	"bufio"
	"fmt"
	"io"
	"os"
	"sync"

	"github.com/pkg/errors"
	"github.com/shenwei356/bio/taxdump"
	"github.com/shenwei356/unik/v5"

	"github.com/spf13/cobra"
	"github.com/twotwotwo/sorts/sortutil"
)

var diffCmd = &cobra.Command{
	Use:   "diff",
	Short: "Set difference of multiple binary files",
	Long: `Set difference of multiple binary files

Attentions:
  0. The first file should be sorted.
  1. The 'canonical/scaled/hashed' flags of all files should be consistent.
  2. By default taxids in the 2nd and later files are ignored.
  3. You can switch on flag -t/--compare-taxid , and input
     files should ALL have or don't have taxid information.
     A same k-mer found but query taxid equals to target taxid,
     or query taxid is ancester of target taxid, this k-mer remains

Tips:
  1. Increasing threads number (-j/--threads) to accelerate computation
     when dealing with lots of files, in cost of more memory occupation.

`,
	Run: func(cmd *cobra.Command, args []string) {
		opt := getOptions(cmd)

		var err error

		if opt.Verbose {
			log.Info("checking input files ...")
		}
		files := getFileListFromArgsAndFile(cmd, args, true, "infile-list", true)
		if opt.Verbose {
			if len(files) == 1 && isStdin(files[0]) {
				log.Info("no files given, reading from stdin")
			} else {
				log.Infof("%d input file(s) given", len(files))
			}
		}

		checkFileSuffix(opt, extDataFile, files...)

		var nfiles = len(files)

		outFile := getFlagString(cmd, "out-prefix")
		sortKmers := getFlagBool(cmd, "sort")
		compareTaxid := getFlagBool(cmd, "compare-taxid")

		threads := opt.NumCPUs

		mc := make([]CodeTaxid, 0, mapInitSize)

		var infh *bufio.Reader
		var r *os.File
		var reader0 *unik.Reader
		var code uint64
		var taxid uint32
		var k int = -1
		var canonical bool
		var hashed bool
		var hasTaxid bool
		var ok bool

		var taxondb *taxdump.Taxonomy

		// -----------------------------------------------------------------------

		file := files[0]
		if opt.Verbose {
			log.Infof("processing file (%d/%d): %s", 1, nfiles, file)
		}

		// read firstFile

		infh, r, _, err = inStream(file)
		checkError(err)

		reader, err := unik.NewReader(infh)
		checkError(errors.Wrap(err, file))

		if !reader.IsSorted() { // query is sorted
			checkError(fmt.Errorf("the first file should be sorted"))
		}

		reader0 = reader
		k = reader.K
		canonical = reader.IsCanonical()
		hashed = reader.IsHashed()
		hasTaxid = !opt.IgnoreTaxid && reader.HasTaxidInfo()
		if compareTaxid {
			if hasTaxid {
				if opt.Verbose {
					log.Infof("taxids found in file: %s", file)
				}
				taxondb = loadTaxonomy(opt, false)
			} else {
				log.Warningf("not taxids found in file: %s, flag -t/--compare-taxid ignored", file)
			}
		}

		var n0 int
		for {
			code, taxid, err = reader.ReadCodeWithTaxid()
			if err != nil {
				if err == io.EOF {
					break
				}
				checkError(errors.Wrap(err, file))
			}

			mc = append(mc, CodeTaxid{Code: code, Taxid: taxid})
		}
		n0 = len(mc)

		r.Close()

		if opt.Verbose {
			log.Infof("%d k-mers loaded", n0)
		}

		if n0 == 0 {
			if opt.Verbose {
				log.Infof("exporting k-mers")
			}

			if !isStdout(outFile) {
				outFile += extDataFile
			}
			outfh, gw, w, err := outStream(outFile, opt.Compress, opt.CompressionLevel)
			checkError(err)
			defer func() {
				outfh.Flush()
				if gw != nil {
					gw.Close()
				}
				w.Close()
			}()

			var mode uint32
			if sortKmers {
				mode |= unik.UnikSorted
			} else if opt.Compact && !hashed {
				mode |= unik.UnikCompact
			}
			if canonical {
				mode |= unik.UnikCanonical
			}
			if hasTaxid {
				mode |= unik.UnikIncludeTaxID
			}
			if hashed {
				mode |= unik.UnikHashed
			}

			writer, err := unik.NewWriter(outfh, k, mode)
			checkError(errors.Wrap(err, outFile))
			writer.SetMaxTaxid(maxUint32N(reader.GetTaxidBytesLength())) // follow reader

			writer.Number = 0
			checkError(writer.WriteHeader())
			checkError(writer.Flush())

			if opt.Verbose {
				log.Infof("%d k-mers saved to %s", 0, outFile)
			}
			return
		}

		// -----------------------------------------------------------------------

		if threads > len(files)-1 {
			threads = len(files) - 1
		}
		if threads < 1 {
			threads = 1
		}

		done := make(chan int)

		toStop := make(chan int, threads+2)
		doneDone := make(chan int)
		go func() {
			<-toStop
			close(done)
			doneDone <- 1
		}()

		// ---------------

		type iFile struct {
			i    int
			file string
		}

		chFile := make(chan iFile, threads)
		doneSendFile := make(chan int)

		maps := make(map[int]map[uint64]uint32, threads)

		mapsc := make(map[int][]CodeTaxid, threads)
		mapsc[0] = mc

		if threads > 1 {
			// clone maps
			if opt.Verbose {
				log.Infof("cloning data for parallization")
			}
			var wg sync.WaitGroup
			type iMap struct {
				i  int
				m  map[uint64]uint32
				mc []CodeTaxid
			}
			ch := make(chan iMap, threads)
			doneClone := make(chan int)
			go func() {
				for ptr := range ch {
					mapsc[ptr.i] = ptr.mc
				}
				doneClone <- 1
			}()
			for i := 1; i < threads; i++ {
				wg.Add(1)
				go func(i int) {
					mc1 := make([]CodeTaxid, len(mc))
					copy(mc1, mc)
					ch <- iMap{i: i, mc: mc1}
					wg.Done()
				}(i)
			}
			wg.Wait()
			close(ch)
			<-doneClone
			if opt.Verbose {
				log.Infof("done cloning data")
			}
		}

		// -----------------------------------------------------------------------
		if opt.Verbose {
			log.Infof("%d workers in position", threads)
		}

		hasDiff := true
		var wgWorkers sync.WaitGroup
		for i := 0; i < threads; i++ { // workers
			wgWorkers.Add(1)

			go func(i int) {
				defer func() {
					if opt.Verbose {
						log.Infof("worker %02d: finished with %d k-mers", i, len(maps[i]))
					}
					wgWorkers.Done()
				}()

				if opt.Verbose {
					log.Infof("worker %02d: started", i)
				}

				var code uint64
				var qtaxid, taxid uint32
				var ifile iFile
				var file string
				var infh *bufio.Reader
				var r *os.File
				var reader *unik.Reader
				var ok bool
				var sorted bool
				var m1 map[uint64]uint32
				mc1 := mapsc[i]
				for {
					ifile, ok = <-chFile
					if !ok {
						return
					}
					file = ifile.file

					select {
					case <-done:
						return
					default:
					}

					if opt.Verbose {
						log.Infof("worker %02d: starting processing file (%d/%d): %s", i, ifile.i+1, nfiles, file)
					}

					infh, r, _, err = inStream(file)
					checkError(err)

					reader, err = unik.NewReader(infh)
					checkError(errors.Wrap(err, file))

					checkCompatibility(reader0, reader, file)
					if compareTaxid && reader.HasTaxidInfo() != hasTaxid {
						if reader.HasTaxidInfo() {
							checkError(fmt.Errorf(`taxid information not found in previous files, but found in this: %s`, file))
						} else {
							checkError(fmt.Errorf(`taxid information found in previous files, but missing in this: %s`, file))
						}
					}

					// file is sorted, so we can skip codes that are small than minCode
					sorted = reader.IsSorted()

					if !sorted {
						if m1 == nil { // clone mc for this worker
							m1 = make(map[uint64]uint32, len(mc))
							for _, ct := range mc {
								m1[ct.Code] = ct.Taxid
							}
							maps[i] = m1
						}

						for {
							code, taxid, err = reader.ReadCodeWithTaxid()
							if err != nil {
								if err == io.EOF {
									break
								}
								checkError(errors.Wrap(err, file))
							}

							// delete seen kmer
							if qtaxid, ok = m1[code]; ok { // slowest part
								if compareTaxid && (qtaxid == taxid ||
									taxondb.LCA(taxid, qtaxid) == qtaxid) {
									continue
								}
								delete(m1, code)
							}
						}

						r.Close()

						if opt.Verbose {
							log.Infof("worker %02d: finished processing file (%d/%d): %s, %d k-mers remain", i, ifile.i+1, nfiles, file, len(m1))
						}
						if len(m1) == 0 {
							hasDiff = false
							toStop <- 1
							return
						}
					} else {
						mc2 := make([]CodeTaxid, 0, len(mc1))
						var qCode, code uint64
						var qtaxid, taxid uint32
						ii := 0

						qCode = mc1[ii].Code
						qtaxid = mc1[ii].Taxid
						code, taxid, err = reader.ReadCodeWithTaxid()
						if err != nil {
							if err == io.EOF {
								break
							}
							checkError(errors.Wrap(err, file))
						}

						for {
							if qCode < code {
								mc2 = append(mc2, mc1[ii])

								ii++
								if ii >= len(mc1) {
									break
								}
								qCode = mc1[ii].Code
								qtaxid = mc1[ii].Taxid
							} else if qCode == code {
								if compareTaxid && (qtaxid == taxid || // keep k-mer with same taxid
									taxondb.LCA(taxid, qtaxid) == qtaxid) { // keep k-mer which is son of query
									mc2 = append(mc2, mc1[ii])
								}

								ii++
								if ii >= len(mc1) {
									break
								}
								qCode = mc1[ii].Code
								qtaxid = mc1[ii].Taxid

								code, taxid, err = reader.ReadCodeWithTaxid()
								if err != nil {
									if err == io.EOF {
										break
									}
									checkError(errors.Wrap(err, file))
								}
							} else {
								code, taxid, err = reader.ReadCodeWithTaxid()
								if err != nil {
									if err == io.EOF {
										break
									}
									checkError(errors.Wrap(err, file))
								}
							}
						}
						mc2 = append(mc2, mc1[ii:]...)

						r.Close()

						mc1 = mc2
						if opt.Verbose {
							log.Infof("worker %02d: finished processing file (%d/%d): %s, %d k-mers remain", i, ifile.i+1, nfiles, file, len(mc1))
						}
						if len(mc1) == 0 {
							hasDiff = false
							toStop <- 1
							return
						}

						m1 = make(map[uint64]uint32, len(mc1))
						for _, ct := range mc1 {
							m1[ct.Code] = ct.Taxid
						}
						maps[i] = m1
					}

				}
			}(i)
		}

		// send file
		go func() {
		SENDFILE:
			for i, file := range files[1:] {
				if file == files[0] {
					continue
				}
				select {
				case <-done:
					break SENDFILE
				default:
				}

				chFile <- iFile{i + 1, file}
			}
			close(chFile)

			doneSendFile <- 1
		}()

		<-doneSendFile
		wgWorkers.Wait()
		toStop <- 1
		<-doneDone

		var m0 map[uint64]uint32
		if !hasDiff {
			if opt.Verbose {
				log.Infof("no set difference found")
			}
			// return
		} else {
			if opt.Verbose {
				log.Infof("merging results from workers")
			}
			var code uint64
			for _, m := range maps {
				if len(m) == 0 {
					m0 = m
					break
				}

				if m0 == nil {
					m0 = m
					continue
				}
				for code = range m0 {
					if _, ok = m[code]; !ok { // it's already been deleted in other m
						delete(m0, code) // so it should be deleted
					}
				}

				if len(m0) == 0 {
					break
				}
			}

			if len(m0) == 0 {
				if opt.Verbose {
					log.Warningf("no set difference found")
				}
				// return
			}
		}

		// -----------------------------------------------------------------------

		// output

		if opt.Verbose {
			log.Infof("exporting k-mers")
		}

		if !isStdout(outFile) {
			outFile += extDataFile
		}
		outfh, gw, w, err := outStream(outFile, opt.Compress, opt.CompressionLevel)
		checkError(err)
		defer func() {
			outfh.Flush()
			if gw != nil {
				gw.Close()
			}
			w.Close()
		}()

		var mode uint32
		if sortKmers {
			mode |= unik.UnikSorted
		} else if opt.Compact && !hashed {
			mode |= unik.UnikCompact
		}
		if canonical {
			mode |= unik.UnikCanonical
		}
		if hasTaxid {
			mode |= unik.UnikIncludeTaxID
		}
		if hashed {
			mode |= unik.UnikHashed
		}

		writer, err := unik.NewWriter(outfh, k, mode)
		checkError(errors.Wrap(err, outFile))
		writer.SetMaxTaxid(opt.MaxTaxid)

		if sortKmers {
			writer.Number = uint64(len(m0))
		}

		if len(m0) == 0 {
			writer.Number = 0
			checkError(writer.WriteHeader())
		} else {
			if sortKmers {
				codes := make([]uint64, len(m0))

				i := 0
				for code = range m0 {
					codes[i] = code
					i++
				}

				if opt.Verbose {
					log.Infof("sorting %d k-mers", len(codes))
				}
				// sort.Sort(unikmer.CodeSlice(codes))
				sortutil.Uint64s(codes)
				if opt.Verbose {
					log.Infof("done sorting")
				}

				for _, code = range codes {
					writer.WriteCodeWithTaxid(code, m0[code])
				}
			} else {
				for code, taxid = range m0 {
					writer.WriteCodeWithTaxid(code, taxid)
				}
			}
		}
		checkError(writer.Flush())
		if opt.Verbose {
			log.Infof("%d k-mers saved to %s", len(m0), outFile)
		}
	},
}

func init() {
	RootCmd.AddCommand(diffCmd)

	diffCmd.Flags().StringP("out-prefix", "o", "-", `out file prefix ("-" for stdout)`)
	diffCmd.Flags().BoolP("sort", "s", false, helpSort)
	diffCmd.Flags().BoolP("compare-taxid", "t", false, `take taxid into consideration. type unikmer "diff -h" for detail`)
}
