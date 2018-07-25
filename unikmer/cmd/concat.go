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

	"github.com/shenwei356/unikmer"
	"github.com/shenwei356/xopen"
	"github.com/spf13/cobra"
)

// concatCmd represents
var concatCmd = &cobra.Command{
	Use:   "concat",
	Short: "concatenate multiple binary files",
	Long: `concatenate multiple binary files

`,
	Run: func(cmd *cobra.Command, args []string) {
		opt := getOptions(cmd)
		runtime.GOMAXPROCS(opt.NumCPUs)
		files := getFileList(args)

		outFile := getFlagString(cmd, "out-prefix")

		var err error

		if !isStdout(outFile) {
			outFile += extDataFile
		}
		outfh, err := xopen.WopenGzip(outFile)
		checkError(err)
		defer outfh.Close()

		var writer *unikmer.Writer

		var infh *xopen.Reader
		var reader *unikmer.Reader
		var kcode unikmer.KmerCode
		var k int = -1
		var flag int
		var nfiles = len(files)
		for i, file := range files {
			if !isStdin(file) && !strings.HasSuffix(file, extDataFile) {
				checkError(fmt.Errorf("input should be stdin or %s file", extDataFile))
			}

			if opt.Verbose {
				log.Infof("process file (%d/%d): %s", i+1, nfiles, file)
			}

			flag = func() int {
				infh, err = xopen.Ropen(file)
				checkError(err)
				defer infh.Close()

				reader, err = unikmer.NewReader(infh)
				checkError(err)

				if k == -1 {
					k = reader.K
					writer = unikmer.NewWriter(outfh, k)
				} else if k != reader.K {
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

					writer.Write(kcode) // not need to check err
				}

				return flagContinue
			}()

			if flag == flagReturn {
				return
			} else if flag == flagBreak {
				break
			}
		}

	},
}

func init() {
	RootCmd.AddCommand(concatCmd)

	concatCmd.Flags().StringP("out-prefix", "o", "-", `out file prefix ("-" for stdout)`)
}
