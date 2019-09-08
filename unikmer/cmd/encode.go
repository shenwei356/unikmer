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
	"fmt"
	"runtime"
	"strings"

	"github.com/shenwei356/breader"
	"github.com/shenwei356/unikmer"
	"github.com/spf13/cobra"
)

// encodeCmd represents
var encodeCmd = &cobra.Command{
	Use:   "encode",
	Short: "encode plain k-mer text to integer",
	Long: `encode plain k-mer text to integer

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

		checkFiles("", files...)

		outFile := getFlagString(cmd, "out-file")
		all := getFlagBool(cmd, "all")
		canonical := getFlagBool(cmd, "canonical")

		outfh, gw, w, err := outStream(outFile, strings.HasSuffix(strings.ToLower(outFile), ".gz"), opt.CompressionLevel)
		checkError(err)
		defer func() {
			outfh.Flush()
			if gw != nil {
				gw.Close()
			}
			w.Close()
		}()

		var k int = -1
		var l int
		var reader *breader.BufferedReader
		var chunk breader.Chunk
		var data interface{}
		var line string
		var kcode unikmer.KmerCode

		for _, file := range files {
			reader, err = breader.NewDefaultBufferedReader(file)
			checkError(err)

			for chunk = range reader.Ch {
				checkError(chunk.Err)
				for _, data = range chunk.Data {
					line = data.(string)
					l = len(line)

					if l == 0 {
						continue
					} else if k == -1 {
						k = l
					} else if l != k {
						checkError(fmt.Errorf("K-mer length mismatch, previous: %d, current: %d. %s", k, l, line))
					}

					kcode, err = unikmer.NewKmerCode([]byte(line))
					if err != nil {
						checkError(fmt.Errorf("fail to encode '%s': %s", line, err))
					}
					if canonical {
						kcode = kcode.Canonical()
					}

					if all {
						outfh.WriteString(fmt.Sprintf("%s\t%s\t%d\t%s\n", line, kcode.String(), kcode.Code, kcode.BitsString()))
					} else {
						outfh.WriteString(fmt.Sprintf("%d\n", kcode.Code))
					}
				}
			}
		}

	},
}

func init() {
	RootCmd.AddCommand(encodeCmd)

	encodeCmd.Flags().StringP("out-file", "o", "-", `out file ("-" for stdout, suffix .gz for gzipped out)`)
	encodeCmd.Flags().BoolP("all", "a", false, `output all data: orginial k-mer, parsed k-mer, encoded integer, encode bits`)
	encodeCmd.Flags().BoolP("canonical", "K", false, "keep the canonical k-mers")
}
