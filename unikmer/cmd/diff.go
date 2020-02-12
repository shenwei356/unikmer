// Copyright © 2018-2020 Wei Shen <shenwei356@gmail.com>
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
	"runtime"
	"sort"
	"sync"

	"github.com/shenwei356/unikmer"
	"github.com/spf13/cobra"
)

// diffCmd represents
var diffCmd = &cobra.Command{
	Use:   "diff",
	Short: "set difference of multiple binary files",
	Long: `set difference of multiple binary files

Attentions:
  1. the 'canonical' flags of all files should be consistent.

Tips:
  1. Increasing threads number (-j/--threads) to accelerate computation,
     in cost of more memory occupation.

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

		checkFileSuffix(extDataFile, files...)

		var nfiles = len(files)
		if nfiles <= 1 {
			checkError(fmt.Errorf("two or more files needed"))
		}

		outFile := getFlagString(cmd, "out-prefix")
		sortKmers := getFlagBool(cmd, "sort")

		threads := opt.NumCPUs

		runtime.GOMAXPROCS(threads)

		m := make(map[uint64]uint32, mapInitSize)

		var infh *bufio.Reader
		var r *os.File
		var reader *unikmer.Reader
		var code uint64
		var taxid uint32
		var k int = -1
		var canonical bool
		var hasTaxid bool
		var ok bool

		// -----------------------------------------------------------------------

		file := files[0]
		if opt.Verbose {
			log.Infof("processing file (%d/%d): %s", 1, nfiles, file)
		}

		// read firstFile

		infh, r, _, err = inStream(file)
		checkError(err)

		reader, err = unikmer.NewReader(infh)
		checkError(err)

		k = reader.K
		canonical = reader.IsCanonical()
		hasTaxid = !opt.IgnoreTaxid && reader.HasTaxidInfo()

		var minCode uint64 = ^uint64(0)
		if reader.IsSorted() { // query is sorted
			once := true
			for {
				code, taxid, err = reader.ReadCodeWithTaxid()
				if err != nil {
					if err == io.EOF {
						break
					}
					checkError(err)
				}

				if once {
					minCode = code
					once = false
				}
				m[code] = taxid
			}
		} else {
			for {
				code, taxid, err = reader.ReadCodeWithTaxid()
				if err != nil {
					if err == io.EOF {
						break
					}
					checkError(err)
				}

				if code < minCode {
					minCode = code
				}
				m[code] = taxid
			}
		}

		r.Close()

		if opt.Verbose {
			log.Infof("%d k-mers loaded", len(m))
			log.Infof("min code: %s (%d)", unikmer.KmerCode{Code: minCode, K: k}, minCode)
		}

		if len(m) == 0 {
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
				mode |= unikmer.UNIK_SORTED
			} else if opt.Compact {
				mode |= unikmer.UNIK_COMPACT
			}
			if canonical {
				mode |= unikmer.UNIK_CANONICAL
			}
			if hasTaxid {
				mode |= unikmer.UNIK_INCLUDETAXID
			}

			writer, err := unikmer.NewWriter(outfh, k, mode)
			checkError(err)

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
		maps[0] = m

		if len(files) > 2 {
			// clone maps
			if opt.Verbose {
				log.Infof("cloning data for parallization")
			}
			var wg sync.WaitGroup
			type iMap struct {
				i int
				m map[uint64]uint32
			}
			ch := make(chan iMap, threads)
			doneClone := make(chan int)
			go func() {
				for ptr := range ch {
					maps[ptr.i] = ptr.m
				}
				doneClone <- 1
			}()
			for i := 1; i < threads; i++ {
				wg.Add(1)
				go func(i int) {
					m1 := make(map[uint64]uint32, len(m))
					for k, t := range m {
						m1[k] = t
					}
					ch <- iMap{i: i, m: m1}
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
		log.Infof("%d workers in position", threads)

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
				var ifile iFile
				var file string
				var infh *bufio.Reader
				var r *os.File
				var reader *unikmer.Reader
				var ok bool
				var sorted bool
				var nSkip int
				var checkSkip bool
				m1 := maps[i]
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

					reader, err = unikmer.NewReader(infh)
					checkError(err)

					if k != reader.K {
						checkError(fmt.Errorf("K (%d) of binary file '%s' not equal to previous K (%d)", reader.K, file, k))
					}
					if reader.IsCanonical() != canonical {
						checkError(fmt.Errorf(`'canonical' flags not consistent, please check with "unikmer stats"`))
					}
					// if !opt.IgnoreTaxid && reader.HasTaxidInfo() != hasTaxid {
					// 	checkError(fmt.Errorf(`taxid information found in some files but missing in others, please check with "unikmer stats"`))
					// }

					// file is sorted, so we can skip codes that are small than minCode
					sorted = reader.IsSorted()
					nSkip = 0
					checkSkip = sorted
					for {
						code, _, err = reader.ReadCodeWithTaxid()
						if err != nil {
							if err == io.EOF {
								break
							}
							checkError(err)
						}

						if checkSkip {
							if code < minCode {
								nSkip++
								continue
							} else {
								if opt.Verbose {
									log.Infof("worker %02d: started processing file (%d/%d): leading %d k-mers skipped for comparison", i, ifile.i+1, nfiles, nSkip)
								}
								checkSkip = false
							}
						}

						// delete seen kmer
						if _, ok = m1[code]; ok { // slowest part
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
			mode |= unikmer.UNIK_SORTED
		} else if opt.Compact {
			mode |= unikmer.UNIK_COMPACT
		}
		if canonical {
			mode |= unikmer.UNIK_CANONICAL
		}
		if hasTaxid {
			mode |= unikmer.UNIK_INCLUDETAXID
		}

		writer, err := unikmer.NewWriter(outfh, k, mode)
		checkError(err)

		if sortKmers {
			writer.Number = int64(len(m0))
		}

		if len(m0) == 0 {
			writer.Number = 0
			checkError(writer.WriteHeader())
		} else {
			if sortKmers {
				codes := make([]unikmer.CodeTaxid, len(m0))
				i := 0
				for code, taxid := range m0 {
					codes[i] = unikmer.CodeTaxid{Code: code, Taxid: taxid}
					i++
				}
				if opt.Verbose {
					log.Infof("sorting %d k-mers", len(codes))
				}
				sort.Sort(unikmer.CodeTaxidSlice(codes))
				if opt.Verbose {
					log.Infof("done sorting")
				}
				for _, codeT := range codes {
					writer.WriteCodeWithTaxid(codeT.Code, codeT.Taxid)
				}
			} else {
				for code, taxid := range m0 {
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
}
