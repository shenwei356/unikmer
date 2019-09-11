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

	"github.com/shenwei356/util/pathutil"

	"github.com/shenwei356/unikmer"
	"github.com/spf13/cobra"
)

// splitCmd represents
var splitCmd = &cobra.Command{
	Use:   "split",
	Short: "split k-mers into sorted chunk files",
	Long: `split k-mers into sorted chunk files

`,
	Run: func(cmd *cobra.Command, args []string) {
		opt := getOptions(cmd)
		runtime.GOMAXPROCS(opt.NumCPUs)

		var err error

		var files []string
		infileList := getFlagString(cmd, "infile-list")
		if infileList != "" {
			files, err = getListFromFile(infileList)
			checkError(err)
		} else {
			files = getFileList(args)
		}

		checkFiles(extDataFile, files...)

		outDir := getFlagString(cmd, "out-dir")

		maxMem, err := ParseByteSize(getFlagString(cmd, "chunk-size"))
		if err != nil {
			checkError(fmt.Errorf("parsing byte size: %s", err))
		}
		if maxMem == 0 {
			checkError(fmt.Errorf("non-zero chunk size needed"))
		}
		maxElem := maxMem >> 3 // uint64 == 8 bytes
		if maxMem > 0 && maxElem < 1 {
			maxElem = 1
		}
		m := make([]uint64, 0, mapInitSize)

		existed, err := pathutil.DirExists(outDir)
		checkError(err)
		if !existed {
			err = os.MkdirAll(outDir, 0755)
			if err != nil {
				checkError(fmt.Errorf("fail to create out directory: %s", outDir))
			}
		}

		var infh *bufio.Reader
		var r *os.File
		var reader *unikmer.Reader
		var kcode unikmer.KmerCode
		var k int = -1
		var canonical bool
		var mode uint32
		var firstFile = true
		var flag int
		var nfiles = len(files)

		var iTmpFile int
		limitMem := maxElem > 0

		var wg sync.WaitGroup
		tokens := make(chan int, opt.NumCPUs)

		// just for counting total k-mers
		chN := make(chan int64, opt.NumCPUs)
		var N int64 = 0
		done := make(chan int)
		go func() {
			for n := range chN {
				N += n
			}
			done <- 1
		}()

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

					if opt.Compact {
						mode |= unikmer.UNIK_COMPACT
					}
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
					kcode, err = reader.Read()
					if err != nil {
						if err == io.EOF {
							break
						}
						checkError(err)
					}

					m = append(m, kcode.Code)

					if limitMem && len(m) >= maxElem {
						iTmpFile++
						outFile1 := chunFileName(outDir, iTmpFile)

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
								log.Infof("[chunk %d] %d k-mers saved to %s", iTmpFile, len(m), outFile)
							}
							chN <- int64(len(m))
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

		// dump remaining k-mers to file
		if len(m) > 0 {
			iTmpFile++
			outFile1 := chunFileName(outDir, iTmpFile)

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
					log.Infof("[chunk %d] %d k-mers saved to %s", iTmpFile, len(m), outFile)
				}
				chN <- int64(len(m))
			}(m, iTmpFile, outFile1)
		}

		// wait all k-mers being wrote to files
		wg.Wait()
		close(chN)
		<-done

		if opt.Verbose {
			log.Infof("%d chunk files with total %d k-mers saved to dir: %s", iTmpFile, N, outDir)
		}
	},
}

func init() {
	RootCmd.AddCommand(splitCmd)

	splitCmd.Flags().StringP("out-dir", "o", "./", `output directory`)
	splitCmd.Flags().StringP("chunk-size", "m", "", `split input into chunks of N bytes, supports K/M/G suffix, type "unikmer split -h" for detail`)
}
