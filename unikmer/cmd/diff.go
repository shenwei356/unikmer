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
	"fmt"
	"io"
	"runtime"
	"strings"
	"sync"

	"github.com/shenwei356/unikmer"
	"github.com/shenwei356/xopen"
	"github.com/spf13/cobra"
)

// diffCmd represents
var diffCmd = &cobra.Command{
	Use:   "diff",
	Short: "set difference of multiple binary files",
	Long: `set difference of multiple binary files

`,
	Run: func(cmd *cobra.Command, args []string) {
		opt := getOptions(cmd)
		files := getFileList(args)

		outFile := getFlagString(cmd, "out-prefix")
		// checkInterval := getFlagPositiveInt(cmd, "check-interval")
		threads := opt.NumCPUs

		runtime.GOMAXPROCS(threads)

		var err error

		m := make(map[uint64]bool, mapInitSize)

		var infh *xopen.Reader
		var reader *unikmer.Reader
		var kcode unikmer.KmerCode
		var k int = -1
		var ok bool
		var nfiles = len(files)

		// -----------------------------------------------------------------------

		file := files[0]
		if !isStdin(file) && !strings.HasSuffix(file, extDataFile) {
			checkError(fmt.Errorf("input should be stdin or %s file", extDataFile))
		}
		if opt.Verbose {
			log.Infof("process file (%d/%d): %s", 1, nfiles, file)
		}

		// only one file given
		if len(files) == 1 {
			func() {
				infh, err = xopen.Ropen(file)
				checkError(err)
				defer infh.Close()

				if !isStdout(outFile) {
					outFile += extDataFile
				}

				var outfh *xopen.Writer
				outfh, err = xopen.WopenGzip(outFile)
				checkError(err)
				defer outfh.Close()

				writer := unikmer.NewWriter(outfh, k)

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

		infh, err = xopen.Ropen(file)
		checkError(err)

		reader, err = unikmer.NewReader(infh)
		checkError(err)

		k = reader.K

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

		infh.Close()

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
				var infh *xopen.Reader
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

					infh, err = xopen.Ropen(file)
					checkError(err)

					reader, err = unikmer.NewReader(infh)
					checkError(err)

					if k != reader.K {
						checkError(fmt.Errorf("K (%d) of binary file '%s' not equal to previous K (%d)", reader.K, file, k))
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

					infh.Close()

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
				if !isStdin(file) && !strings.HasSuffix(file, extDataFile) {
					checkError(fmt.Errorf("input should be stdin or %s file", extDataFile))
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
		outfh, err := xopen.WopenGzip(outFile)
		checkError(err)
		defer outfh.Close()

		writer := unikmer.NewWriter(outfh, k)

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
	diffCmd.Flags().IntP("check-interval", "i", 5, `check kmers every N files, N > 1 could save some time`)
}
