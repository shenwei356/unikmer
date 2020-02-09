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

	"github.com/shenwei356/unikmer"
	"github.com/spf13/cobra"
)

// sampleCmd represents
var sampleCmd = &cobra.Command{
	Use:   "sample",
	Short: "sample k-mers from binary files",
	Long: `sample k-mers from binary files.

The sampling type is fixed sampling.

Attentions:
  1. the 'canonical' flags of all files should be consistent.

`,
	Run: func(cmd *cobra.Command, args []string) {
		opt := getOptions(cmd)
		runtime.GOMAXPROCS(opt.NumCPUs)

		var err error

		files := getFileListFromArgsAndFile(cmd, args, true, "infile-list", true)

		checkFileSuffix(extDataFile, files...)

		outFile := getFlagString(cmd, "out-prefix")

		sampling := true
		start := getFlagPositiveInt(cmd, "start")
		window := getFlagPositiveInt(cmd, "window")
		if start == 1 && window == 1 {
			sampling = false
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
		var flag int
		var nfiles = len(files)
		var n uint64
		var j int
		for i, file := range files {
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
					writer, err = unikmer.NewWriter(outfh, k, reader.Flag)
					checkError(err)
				} else if k != reader.K {
					checkError(fmt.Errorf("K (%d) of binary file '%s' not equal to previous K (%d)", reader.K, file, k))
				} else if (reader.Flag&unikmer.UNIK_CANONICAL > 0) != canonical {
					checkError(fmt.Errorf(`'canonical' flags not consistent, please check with "unikmer stats"`))
				}

				if sampling {
					j = 0
					for {
						kcode, err = reader.Read()
						if err != nil {
							if err == io.EOF {
								break
							}
							checkError(err)
						}

						j++
						if (j-start)%window == 0 || j == start {
							n++
							writer.WriteCode(kcode.Code) // not need to check err
						}
					}

				} else {
					for {
						kcode, err = reader.Read()
						if err != nil {
							if err == io.EOF {
								break
							}
							checkError(err)
						}

						n++
						writer.WriteCode(kcode.Code) // not need to check err
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

		checkError(writer.Flush())
		if opt.Verbose {
			log.Infof("%d k-mers saved to %s", n, outFile)
		}
	},
}

func init() {
	RootCmd.AddCommand(sampleCmd)

	sampleCmd.Flags().IntP("start", "s", 1, `start location`)
	sampleCmd.Flags().IntP("window", "w", 1, `window size`)

	sampleCmd.Flags().StringP("out-prefix", "o", "-", `out file prefix ("-" for stdout)`)
}
