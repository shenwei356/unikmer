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
  1. Increasing threads number (-t/--threads) to accelerate computation,
     in cost of more memory occupation.

`,
	Run: func(cmd *cobra.Command, args []string) {
		opt := getOptions(cmd)
		files := getFileList(args)

		checkFiles(files)

		outFile := getFlagString(cmd, "out-prefix")
		threads := opt.NumCPUs

		runtime.GOMAXPROCS(threads)

		var err error

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

		toStop := make(chan int, 1)
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
					select {
					case <-done:
						return
					default:
					}

					ifile, ok = <-chFile
					if !ok {
						return
					}
					file = ifile.file

					if opt.Verbose {
						log.Infof("(worker %d) process file (%d/%d): %s", i, ifile.i+1, nfiles, file)
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
			for i, file := range files[1:] {
				if file == files[0] {
					continue
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

		if !hasDiff {
			if opt.Verbose {
				log.Infof("no set difference found")
			}
			return
		}

		var m0 map[uint64]bool
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
				log.Infof("no set difference found")
			}
			return
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
		writer, err := unikmer.NewWriter(outfh, k, mode)
		checkError(err)

		for code := range m0 {
			writer.Write(unikmer.KmerCode{Code: code, K: k})
		}
		if opt.Verbose {
			log.Infof("%d Kmers saved", len(m0))
		}
	},
}

func init() {
	RootCmd.AddCommand(diffCmd)

	diffCmd.Flags().StringP("out-prefix", "o", "-", `out file prefix ("-" for stdout)`)
}
