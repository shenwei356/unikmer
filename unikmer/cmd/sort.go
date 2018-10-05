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

	"github.com/shenwei356/unikmer"
	"github.com/spf13/cobra"
)

// sortCmd represents
var sortCmd = &cobra.Command{
	Use:   "sort",
	Short: "sort binary files",
	Long: `sort binary files
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

		checkFiles(files)

		outFile := getFlagString(cmd, "out-prefix")
		unique := getFlagBool(cmd, "unique")

		m := make([]uint64, 0, mapInitSize)

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

		var writer *unikmer.Writer

		var infh *bufio.Reader
		var r *os.File
		var reader *unikmer.Reader
		var kcode unikmer.KmerCode
		var k int = -1
		var canonical bool
		var firstFile = true
		var flag int
		var nfiles = len(files)
		for i, file := range files {
			if !firstFile && file == files[0] {
				continue
			}

			if opt.Verbose {
				log.Infof("process file (%d/%d): %s", i+1, nfiles, file)
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

					var mode uint32
					if opt.Compact {
						mode |= unikmer.UNIK_COMPACT
					}
					if canonical {
						mode |= unikmer.UNIK_CANONICAL
					}
					mode |= unikmer.UNIK_SORTED
					writer, err = unikmer.NewWriter(outfh, k, mode)
					checkError(err)
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

				}

				return flagContinue
			}()

			if flag == flagReturn {
				return
			} else if flag == flagBreak {
				break
			}
		}

		sort.Sort(unikmer.CodeSlice(m))

		var n int
		if unique {
			var last uint64 = ^uint64(0)
			for _, code := range m {
				if code == last {
					continue
				}
				last = code
				n++
				writer.Write(unikmer.KmerCode{code, k})
			}
		} else {
			for _, code := range m {
				writer.Write(unikmer.KmerCode{code, k})
			}
			n = len(m)
		}

		if opt.Verbose {
			log.Infof("%d Kmers saved", n)
		}
	},
}

func init() {
	RootCmd.AddCommand(sortCmd)

	sortCmd.Flags().StringP("out-prefix", "o", "-", `out file prefix ("-" for stdout)`)
	sortCmd.Flags().BoolP("unique", "u", false, `remove duplicated Kmers`)
}
