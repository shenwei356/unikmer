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

// unionCmd represents
var unionCmd = &cobra.Command{
	Use:     "union",
	Aliases: []string{"uniq"},
	Short:   "union of multiple binary files",
	Long: `union of multiple binary files

Attentions:
  1. the 'canonical' flags of all files should be consistent.

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

		outFile := getFlagString(cmd, "out-prefix")
		sortKmers := getFlagBool(cmd, "sort")
		repeated := getFlagBool(cmd, "repeated")

		var m map[uint64]struct{}
		var mb map[uint64]bool
		if repeated {
			mb = make(map[uint64]bool, mapInitSize)
		} else {
			m = make(map[uint64]struct{}, mapInitSize)
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

		var writer *unikmer.Writer

		var infh *bufio.Reader
		var r *os.File
		var reader *unikmer.Reader
		var kcode unikmer.KmerCode
		var k int = -1
		var canonical bool
		var firstFile = true
		var existed, ok bool
		var n int64
		var flag int
		var nfiles = len(files)
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

					if !sortKmers {
						var mode uint32
						if opt.Compact {
							mode |= unikmer.UNIK_COMPACT
						}
						if canonical {
							mode |= unikmer.UNIK_CANONICAL
						}
						writer, err = unikmer.NewWriter(outfh, k, mode)
						checkError(err)
					}
				} else if k != reader.K {
					checkError(fmt.Errorf("K (%d) of binary file '%s' not equal to previous K (%d)", reader.K, file, k))
				} else if (reader.Flag&unikmer.UNIK_CANONICAL > 0) != canonical {
					checkError(fmt.Errorf(`'canonical' flags not consistent, please check with "unikmer stats"`))
				}

				if repeated {
					for {
						kcode, err = reader.Read()
						if err != nil {
							if err == io.EOF {
								break
							}
							checkError(err)
						}

						// new kmers
						if existed, ok = mb[kcode.Code]; !ok {
							mb[kcode.Code] = false
						} else {
							if !existed {
								mb[kcode.Code] = true // mark repeated
								n++
								if !sortKmers {
									writer.Write(kcode) // not need to check err
								}
							}
						}
					}

					return flagContinue
				}

				for {
					kcode, err = reader.Read()
					if err != nil {
						if err == io.EOF {
							break
						}
						checkError(err)
					}

					// new kmers
					if _, ok = m[kcode.Code]; !ok {
						m[kcode.Code] = struct{}{}
						n++
						if !sortKmers {
							writer.Write(kcode) // not need to check err
						}
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

		if sortKmers {
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

			writer.Number = int64(len(m))

			var codes []uint64
			if repeated {
				codes = make([]uint64, 0, len(m))
			} else {
				codes = make([]uint64, len(m))
			}

			i := 0
			if repeated {
				for code, r := range mb {
					if r {
						codes = append(codes, code)
					}
				}
			} else {
				for code := range m {
					codes[i] = code
					i++
				}
			}
			if opt.Verbose {
				log.Infof("sorting %d k-mers", len(codes))
			}
			sort.Sort(unikmer.CodeSlice(codes))
			if opt.Verbose {
				log.Infof("done sorting")
			}
			for _, code := range codes {
				writer.Write(unikmer.KmerCode{Code: code, K: k})
			}
		}

		checkError(writer.Flush())
		if opt.Verbose {
			log.Infof("%d k-mers saved", n)
		}
	},
}

func init() {
	RootCmd.AddCommand(unionCmd)

	unionCmd.Flags().StringP("out-prefix", "o", "-", `out file prefix ("-" for stdout)`)
	unionCmd.Flags().BoolP("sort", "s", false, helpSort)
	unionCmd.Flags().BoolP("repeated", "d", false, `only print duplicate k-mers`)
}
