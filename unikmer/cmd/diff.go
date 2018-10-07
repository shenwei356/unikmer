// Copyright © 2018 Wei Shen <shenwei356@gmail.com>
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

		var files []string
		infileList := getFlagString(cmd, "infile-list")
		if infileList != "" {
			files, err = getListFromFile(infileList)
			checkError(err)
		} else {
			files = getFileList(args)
		}

		checkFiles(files)

		outFile := getFlagString(cmd, "out-prefix")
		sortKmers := getFlagBool(cmd, "sort")

		threads := opt.NumCPUs

		runtime.GOMAXPROCS(threads)

		m := make(map[uint64]bool, mapInitSize)

		var infh *bufio.Reader
		var r *os.File
		var reader *unikmer.Reader
		var kcode unikmer.KmerCode
		var k int = -1
		var canonical bool
		var ok bool
		var nfiles = len(files)

		// -----------------------------------------------------------------------

		file := files[0]
		if opt.Verbose {
			log.Infof("process file (%d/%d): %s", 1, nfiles, file)
		}

		// only one file given
		if len(files) == 1 {
			func() {
				infh, r, _, err = inStream(file)
				checkError(err)
				defer r.Close()

				reader, err = unikmer.NewReader(infh)
				checkError(err)

				if !isStdout(outFile) {
					outFile += extDataFile
				}
				outfh, gw, w, err := outStream(outFile, opt.Compress)
				checkError(err)
				defer func() {
					outfh.Flush()
					if gw != nil {
						gw.Close()
					}
					w.Close()
				}()

				var mode uint32
				if opt.Compact {
					mode |= unikmer.UNIK_COMPACT
				}
				if canonical {
					mode |= unikmer.UNIK_CANONICAL
				}
				writer, err := unikmer.NewWriter(outfh, reader.K, mode)
				checkError(err)

				m := make(map[uint64]struct{}, mapInitSize)
				for {
					kcode, err = reader.Read()
					if err != nil {
						if err == io.EOF {
							break
						}
						checkError(err)
					}

					if _, ok = m[kcode.Code]; !ok {
						m[kcode.Code] = struct{}{}
						writer.Write(kcode) // not need to check er
					}
				}

				checkError(writer.Flush())
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

		for {
			kcode, err = reader.Read()
			if err != nil {
				if err == io.EOF {
					break
				}
				checkError(err)
			}

			m[kcode.Code] = false
		}

		r.Close()

		if opt.Verbose {
			log.Infof("%d Kmers loaded", len(m))
		}
		// -----------------------------------------------------------------------

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

		maps := make(map[int]map[uint64]bool, threads)
		maps[0] = m

		// clone maps
		if opt.Verbose {
			log.Infof("clone data for parallization")
		}
		var wg sync.WaitGroup
		for i := 1; i < opt.NumCPUs; i++ {
			wg.Add(1)
			go func(i int) {
				m1 := make(map[uint64]bool, len(m))
				for k := range m {
					m1[k] = false
				}
				wg.Done()
				maps[i] = m1
			}(i)
		}
		wg.Wait()
		if opt.Verbose {
			log.Infof("done cloning data")
		}

		// -----------------------------------------------------------------------
		hasDiff := true
		var wgWorkers sync.WaitGroup
		for i := 0; i < opt.NumCPUs; i++ { // workers
			wgWorkers.Add(1)

			go func(i int) {
				defer func() {
					if opt.Verbose {
						log.Infof("worker %d finished with %d Kmers", i, len(maps[i]))
					}
					wgWorkers.Done()
				}()

				if opt.Verbose {
					log.Infof("worker %d started", i)
				}

				var code uint64
				var ifile iFile
				var file string
				var infh *bufio.Reader
				var r *os.File
				var reader *unikmer.Reader
				var kcode unikmer.KmerCode
				var ok, mark bool
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
						log.Infof("(worker %d) process file (%d/%d): %s", i, ifile.i+1, nfiles, file)
					}

					infh, r, _, err = inStream(file)
					checkError(err)

					// if sampling {
					// 	reader, err = unikmer.NewSamplingReader(infh, start, window)
					// } else {
					// 	reader, err = unikmer.NewReader(infh)
					// }
					reader, err = unikmer.NewReader(infh)
					checkError(err)

					if k != reader.K {
						checkError(fmt.Errorf("K (%d) of binary file '%s' not equal to previous K (%d)", reader.K, file, k))
					}

					if (reader.Flag&unikmer.UNIK_CANONICAL > 0) != canonical {
						checkError(fmt.Errorf(`'canonical' flags not consistent, please check with "unikmer stats"`))
					}

					for {
						kcode, err = reader.Read()
						if err != nil {
							if err == io.EOF {
								break
							}
							checkError(err)
						}

						// mark seen kmer
						if _, ok = m1[kcode.Code]; ok {
							m1[kcode.Code] = true
						}
					}

					r.Close()

					// remove seen kmers
					for code, mark = range m1 {
						if mark {
							delete(m1, code)
						}
					}

					if opt.Verbose {
						log.Infof("(worker %d) %d Kmers remain", i, len(m1))
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

		var m0 map[uint64]bool
		if !hasDiff {
			if opt.Verbose {
				log.Infof("no set difference found")
			}
			// return
		} else {
			var code uint64
			for _, m := range maps {
				if len(m) == 0 {
					continue
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
			log.Infof("export Kmers")
		}

		if !isStdout(outFile) {
			outFile += extDataFile
		}
		outfh, gw, w, err := outStream(outFile, opt.Compress)
		checkError(err)
		defer func() {
			outfh.Flush()
			if gw != nil {
				gw.Close()
			}
			w.Close()
		}()

		var mode uint32
		if opt.Compact {
			mode |= unikmer.UNIK_COMPACT
		}
		if canonical {
			mode |= unikmer.UNIK_CANONICAL
		}
		if sortKmers {
			mode |= unikmer.UNIK_SORTED
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
					log.Infof("sort %d Kmers", len(codes))
				}
				sort.Sort(unikmer.CodeSlice(codes))
				if opt.Verbose {
					log.Infof("done sorting")
				}
				for _, code := range codes {
					writer.Write(unikmer.KmerCode{Code: code, K: k})
				}
			} else {
				for code := range m0 {
					writer.Write(unikmer.KmerCode{Code: code, K: k})
				}
			}
		}
		checkError(writer.Flush())
		if opt.Verbose {
			log.Infof("%d Kmers saved", len(m0))
		}
	},
}

func init() {
	RootCmd.AddCommand(diffCmd)

	diffCmd.Flags().StringP("out-prefix", "o", "-", `out file prefix ("-" for stdout)`)
	diffCmd.Flags().BoolP("sort", "s", false, "sort Kmers, this reduces file size, you can even disable gzip compression by flag -C/--no-compress")
}
