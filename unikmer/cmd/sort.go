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
	"path/filepath"
	"runtime"
	"sort"
	"sync"

	"github.com/shenwei356/unikmer"
	"github.com/spf13/cobra"
)

// sortCmd represents
var sortCmd = &cobra.Command{
	Use:   "sort",
	Short: "sort k-mers in binary files to reduce file size",
	Long: `sort k-mers in binary files to reduce file size

Notes:
  1. When sorting from large number of files, this command is equivalent to
     'unikmer split' + 'unikmer merge'.

Tips:
  1. You can use '-m/--chunk-size' to limit memory usage, though which is
     not precise and actually RSS is higher than '-m' * '-j'.
  2. Increasing value of -j/--threads can accelerates splitting stage,
     in cost of more memory occupation.

`,
	Run: func(cmd *cobra.Command, args []string) {
		opt := getOptions(cmd)
		runtime.GOMAXPROCS(opt.NumCPUs)

		var err error

		infileList := getFlagString(cmd, "infile-list")

		files := getFileList(args, true)
		if infileList != "" {
			_files, err := getListFromFile(infileList, true)
			checkError(err)
			files = append(files, _files...)
		}

		checkFileSuffix(extDataFile, files...)

		outFile0 := getFlagString(cmd, "out-prefix")
		unique := getFlagBool(cmd, "unique")
		repeated := getFlagBool(cmd, "repeated")
		tmpDir := getFlagString(cmd, "tmp-dir")
		maxOpenFiles := getFlagPositiveInt(cmd, "max-open-files")

		maxMem, err := ParseByteSize(getFlagString(cmd, "chunk-size"))
		if err != nil {
			checkError(fmt.Errorf("parsing byte size: %s", err))
		}
		maxElem := maxMem >> 3 // uint64 == 8 bytes
		if maxMem > 0 && maxElem < 1 {
			maxElem = 1
		}
		m := make([]uint64, 0, mapInitSize)

		outFile := outFile0
		if !isStdout(outFile) {
			outFile += extDataFile
		}

		var writer *unikmer.Writer

		var infh *bufio.Reader
		var r *os.File
		var reader *unikmer.Reader
		var code uint64
		var k int = -1
		var canonical bool
		var mode uint32
		var firstFile = true
		var flag int
		var nfiles = len(files)

		var tmpFiles []string
		var iTmpFile int
		var hasTmpFile bool
		limitMem := maxElem > 0

		var wg sync.WaitGroup
		tokens := make(chan int, opt.NumCPUs)

		for i, file := range files {
			if !firstFile && file == files[0] {
				continue
			}

			if opt.Verbose {
				log.Infof("processing file (%d/%d): %s", i+1, nfiles, file)
			}

			flag = func() int {
				infh, r, _, err = inStream(file)
				checkError(err)
				defer r.Close()

				reader, err = unikmer.NewReader(infh)
				checkError(err)

				if k == -1 {
					k = reader.K
					canonical = reader.Flag&unikmer.UNIK_CANONICAL > 0

					if canonical {
						mode |= unikmer.UNIK_CANONICAL
					}
					mode |= unikmer.UNIK_SORTED
				} else if k != reader.K {
					checkError(fmt.Errorf("K (%d) of binary file '%s' not equal to previous K (%d)", reader.K, file, k))
				} else if (reader.Flag&unikmer.UNIK_CANONICAL > 0) != canonical {
					checkError(fmt.Errorf(`'canonical' flags not consistent, please check with "unikmer stats"`))
				}

				for {
					code, err = reader.ReadCode()
					if err != nil {
						if err == io.EOF {
							break
						}
						checkError(err)
					}

					m = append(m, code)

					if limitMem && len(m) >= maxElem {
						if !hasTmpFile {
							tmpFiles = make([]string, 0, 10)
							hasTmpFile = true

							if isStdout(outFile0) {
								tmpDir = filepath.Join(tmpDir, "stdout.tmp")
							} else {
								tmpDir = filepath.Join(tmpDir, filepath.Base(outFile0)+".tmp")
							}
							err = os.MkdirAll(tmpDir, 0755)
							if err != nil {
								checkError(fmt.Errorf("fail to create temp directory: %s", tmpDir))
							}
						}
						iTmpFile++
						outFile1 := chunkFileName(tmpDir, iTmpFile)
						tmpFiles = append(tmpFiles, outFile1)

						wg.Add(1)
						tokens <- 1
						go func(m []uint64, iTmpFile int, outFile string) {
							defer func() {
								wg.Done()
								<-tokens
							}()

							if opt.Verbose {
								log.Infof("[chunk %d] sorting %d k-mers", iTmpFile, len(m))
							}
							sort.Sort(unikmer.CodeSlice(m))
							if opt.Verbose {
								log.Infof("[chunk %d] done sorting", iTmpFile)
							}

							dumpCodes2File(m, k, mode, outFile, opt)
							if opt.Verbose {
								log.Infof("[chunk %d] %d k-mers saved to tmp file: %s", iTmpFile, len(m), outFile)
							}
						}(m, iTmpFile, outFile1)

						m = make([]uint64, 0, mapInitSize)
					}

				}

				return flagContinue
			}()

			if flag == flagReturn {
				return
			} else if flag == flagBreak {
				break
			}
		}

		if hasTmpFile {
			// dump remaining k-mers to file
			if len(m) > 0 {
				iTmpFile++
				outFile1 := chunkFileName(tmpDir, iTmpFile)
				tmpFiles = append(tmpFiles, outFile1)

				wg.Add(1)
				tokens <- 1
				go func(m []uint64, iTmpFile int, outFile string) {
					defer func() {
						wg.Done()
						<-tokens
					}()

					if opt.Verbose {
						log.Infof("[chunk %d] sorting %d k-mers", iTmpFile, len(m))
					}
					sort.Sort(unikmer.CodeSlice(m))
					if opt.Verbose {
						log.Infof("[chunk %d] done sorting", iTmpFile)
					}

					dumpCodes2File(m, k, mode, outFile, opt)
					if opt.Verbose {
						log.Infof("[chunk %d] %d k-mers saved to tmp file: %s", iTmpFile, len(m), outFile)
					}
				}(m, iTmpFile, outFile1)
			}

			// wait all k-mers being wrote to files
			wg.Wait()

			// merge sort

			// this implemention was heavily inspired by https://github.com/oxtoacart/emsort/blob/master/emsort.go
			if opt.Verbose {
				log.Infof("merging from %d chunk files", len(tmpFiles))
			}

			files = make([]string, len(tmpFiles))
			copy(files, tmpFiles)
			tmpFiles = make([]string, 0, 10)

			var n int64
			var _files []string
			if len(files) < maxOpenFiles {
				n, _ = mergeChunksFile(opt, files, outFile, k, mode, unique, repeated)
			} else {
				if opt.Verbose {
					log.Infof("two-pass merge performing")
				}

				_files = make([]string, 0, maxOpenFiles)
				for _, file := range files {
					_files = append(_files, file)
					if len(_files) == maxOpenFiles {
						iTmpFile++
						outFile1 := chunkFileName(tmpDir, iTmpFile)

						n, _ := mergeChunksFile(opt, _files, outFile1, k, mode, false, false)
						if opt.Verbose {
							log.Infof("%d k-mers saved to tmp file: %s", n, outFile1)
						}
						tmpFiles = append(tmpFiles, outFile1)
						_files = make([]string, 0, maxOpenFiles)
					}
				}
				if len(_files) > 0 {
					iTmpFile++
					outFile1 := chunkFileName(tmpDir, iTmpFile)

					n, _ := mergeChunksFile(opt, _files, outFile1, k, mode, false, false)
					if opt.Verbose {
						log.Infof("%d k-mers saved to tmp file: %s", n, outFile1)
					}
					tmpFiles = append(tmpFiles, outFile1)
				}

				n, _ = mergeChunksFile(opt, tmpFiles, outFile, k, mode, unique, repeated)
			}
			if opt.Verbose {
				log.Infof("%d k-mers saved to %s", n, outFile)
			}

			// cleanning

			if opt.Verbose {
				log.Infof("removing %d intermediate files", len(tmpFiles)+len(files))
			}
			for _, file := range append(files, tmpFiles...) {
				err := os.Remove(file)
				if err != nil {
					checkError(fmt.Errorf("fail to remove intermediate file: %s", file))
				}
			}
			if opt.Verbose {
				log.Infof("removing tmp dir: %s", tmpDir)
			}
			err = os.Remove(tmpDir)
			if err != nil {
				checkError(fmt.Errorf("fail to remove temp directory: %s", tmpDir))
			}

			return
		}

		// all k-mers are stored in memory

		if opt.Verbose {
			log.Infof("sorting %d k-mers", len(m))
		}
		sort.Sort(unikmer.CodeSlice(m))
		if opt.Verbose {
			log.Infof("done sorting")
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
		writer, err = unikmer.NewWriter(outfh, k, mode)
		checkError(err)

		var n int
		if unique {
			var last uint64 = ^uint64(0)
			for _, code := range m {
				if code == last {
					continue
				}
				last = code
				n++
				writer.WriteCode(code)
			}
		} else if repeated {
			var last uint64 = ^uint64(0)
			var count int
			for _, code := range m {
				if code == last {
					count++
				} else if count > 1 {
					writer.WriteCode(code)
					n++
					last = code
					count = 1
				} else {
					last = code
					count = 1
				}
			}
			if count > 1 {
				writer.WriteCode(code)
				n++
			}
		} else {
			writer.Number = int64(len(m))
			for _, code := range m {
				writer.WriteCode(code)
			}
			n = len(m)
		}

		checkError(writer.Flush())
		if opt.Verbose {
			log.Infof("%d k-mers saved to %s", n, outFile)
		}
	},
}

func init() {
	RootCmd.AddCommand(sortCmd)

	sortCmd.Flags().StringP("out-prefix", "o", "-", `out file prefix ("-" for stdout)`)
	sortCmd.Flags().BoolP("unique", "u", false, `remove duplicated k-mers`)
	sortCmd.Flags().BoolP("repeated", "d", false, `only print duplicate k-mers`)
	sortCmd.Flags().StringP("chunk-size", "m", "", `split input into chunks of N bytes, supports K/M/G suffix, type "unikmer sort -h" for detail`)
	sortCmd.Flags().StringP("tmp-dir", "t", "./", `directory for intermediate files`)
	sortCmd.Flags().IntP("max-open-files", "M", 300, `max number of open files`)
}
