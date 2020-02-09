// Copyright © 2018-2019 Wei Shen <shenwei356@gmail.com>
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

		files := getFileListFromArgsAndFile(cmd, args, true, "infile-list", true)

		checkFileSuffix(extDataFile, files...)

		outFile := getFlagString(cmd, "out-prefix")
		sortKmers := getFlagBool(cmd, "sort")

		threads := opt.NumCPUs

		runtime.GOMAXPROCS(threads)

		m := make(map[uint64]struct{}, mapInitSize)

		var infh *bufio.Reader
		var r *os.File
		var reader *unikmer.Reader
		var code uint64
		var k int = -1
		var canonical bool
		var ok bool
		var nfiles = len(files)

		// -----------------------------------------------------------------------

		file := files[0]
		if opt.Verbose {
			log.Infof("processing file (%d/%d): %s", 1, nfiles, file)
		}

		// only one file given
		if len(files) == 1 {
			func() {
				infh, r, _, err = inStream(file)
				checkError(err)
				defer r.Close()

				reader, err = unikmer.NewReader(infh)
				checkError(err)

				k = reader.K

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

				var writer *unikmer.Writer
				var m2 []uint64

				if sortKmers {
					m2 = make([]uint64, 0, mapInitSize)
				} else {
					var mode uint32
					if sortKmers {
						mode |= unikmer.UNIK_SORTED
					} else if opt.Compact {
						mode |= unikmer.UNIK_COMPACT
					}
					if reader.Flag&unikmer.UNIK_CANONICAL > 0 {
						mode |= unikmer.UNIK_CANONICAL
					}
					writer, err = unikmer.NewWriter(outfh, reader.K, mode)
					checkError(err)
				}

				m := make(map[uint64]struct{}, mapInitSize)
				for {
					code, err = reader.ReadCode()
					if err != nil {
						if err == io.EOF {
							break
						}
						checkError(err)
					}

					if _, ok = m[code]; !ok {
						m[code] = struct{}{}
						if sortKmers {
							m2 = append(m2, code)
						} else {
							writer.WriteCode(code) // not need to check er
						}
					}
				}

				if sortKmers {
					var mode uint32
					if reader.Flag&unikmer.UNIK_CANONICAL > 0 {
						mode |= unikmer.UNIK_CANONICAL
					}
					mode |= unikmer.UNIK_SORTED
					writer, err = unikmer.NewWriter(outfh, reader.K, mode)
					checkError(err)

					writer.Number = int64(len(m2))

					if opt.Verbose {
						log.Infof("sorting %d k-mers", len(m2))
					}
					sort.Sort(unikmer.CodeSlice(m2))
					if opt.Verbose {
						log.Infof("done sorting")
					}

					for _, code = range m2 {
						writer.WriteCode(code)
					}
				}

				checkError(writer.Flush())
				if opt.Verbose {
					log.Infof("%d k-mers saved to %s", len(m), outFile)
				}
			}()

			return
		}

		// -----------------------------------------------------------------------

		// > 1 files

		// read firstFile

		infh, r, _, err = inStream(file)
		checkError(err)

		reader, err = unikmer.NewReader(infh)
		checkError(err)

		k = reader.K
		canonical = reader.Flag&unikmer.UNIK_CANONICAL > 0

		var minCode uint64 = ^uint64(0)
		if reader.Flag&unikmer.UNIK_SORTED > 0 { // query is sorted
			once := true
			for {
				code, err = reader.ReadCode()
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
				m[code] = struct{}{}
			}
		} else {
			for {
				code, err = reader.ReadCode()
				if err != nil {
					if err == io.EOF {
						break
					}
					checkError(err)
				}

				if code < minCode {
					minCode = code
				}
				m[code] = struct{}{}
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

		maps := make(map[int]map[uint64]struct{}, threads)
		maps[0] = m

		if len(files) > 2 {
			// clone maps
			if opt.Verbose {
				log.Infof("cloning data for parallization")
			}
			var wg sync.WaitGroup
			type iMap struct {
				i int
				m map[uint64]struct{}
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
					m1 := make(map[uint64]struct{}, len(m))
					for k := range m {
						m1[k] = struct{}{}
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

					if (reader.Flag&unikmer.UNIK_CANONICAL > 0) != canonical {
						checkError(fmt.Errorf(`'canonical' flags not consistent, please check with "unikmer stats"`))
					}

					// file is sorted, so we can skip codes that are small than minCode
					sorted = reader.Flag&unikmer.UNIK_SORTED > 0
					nSkip = 0
					checkSkip = sorted
					for {
						code, err = reader.ReadCode()
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

		var m0 map[uint64]struct{}
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
				codes := make([]uint64, len(m0))
				i := 0
				for code := range m0 {
					codes[i] = code
					i++
				}
				if opt.Verbose {
					log.Infof("sorting %d k-mers", len(codes))
				}
				sort.Sort(unikmer.CodeSlice(codes))
				if opt.Verbose {
					log.Infof("done sorting")
				}
				for _, code := range codes {
					writer.WriteCode(code)
				}
			} else {
				for code := range m0 {
					writer.WriteCode(code)
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
