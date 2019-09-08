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
	"strconv"
	"strings"

	"github.com/shenwei356/breader"
	"github.com/shenwei356/unikmer"
	"github.com/spf13/cobra"
)

// decodeCmd represents
var decodeCmd = &cobra.Command{
	Use:   "decode",
	Short: "decode encoded integer to k-mer text",
	Long: `decode encoded integer to k-mer text

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
		k := getFlagPositiveInt(cmd, "kmer-len")
		if k > 32 {
			checkError(fmt.Errorf("k > 32 not supported"))
		}

		outfh, gw, w, err := outStream(outFile, strings.HasSuffix(strings.ToLower(outFile), ".gz"), opt.CompressionLevel)
		checkError(err)
		defer func() {
			outfh.Flush()
			if gw != nil {
				gw.Close()
			}
			w.Close()
		}()

		var reader *breader.BufferedReader
		var chunk breader.Chunk
		var data interface{}
		var line string
		var code uint64
		var kmer []byte

		for _, file := range files {
			reader, err = breader.NewDefaultBufferedReader(file)
			checkError(err)

			for chunk = range reader.Ch {
				checkError(chunk.Err)
				for _, data = range chunk.Data {
					line = data.(string)
					if line == "" {
						continue
					}

					code, err = strconv.ParseUint(line, 10, 64)
					if err != nil {
						checkError(fmt.Errorf("encode kmer should be non-negative integer: %s", line))
					}

					if code < 0 {
						checkError(fmt.Errorf("encode kmer should be non-negative integer: %d", code))
					}
					if code > unikmer.MaxCode[k] {
						checkError(fmt.Errorf("encode integer overflows for k=%d (max: %d): %d", k, unikmer.MaxCode[k], code))
					}

					kmer = unikmer.Decode(code, k)
					if err != nil {
						checkError(fmt.Errorf("fail to decode '%s': %s", line, err))
					}

					if all {
						outfh.WriteString(fmt.Sprintf("%d\t%s\n", code, kmer))
					} else {
						outfh.WriteString(fmt.Sprintf("%s\n", kmer))
					}
				}
			}
		}

	},
}

func init() {
	RootCmd.AddCommand(decodeCmd)

	decodeCmd.Flags().StringP("out-file", "o", "-", `out file ("-" for stdout, suffix .gz for gzipped out)`)
	decodeCmd.Flags().IntP("kmer-len", "k", 0, "k-mer length")
	decodeCmd.Flags().BoolP("all", "a", false, `output all data: encoded integer, decoded k-mer`)

}
