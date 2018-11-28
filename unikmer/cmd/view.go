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
	"strings"

	"github.com/shenwei356/unikmer"
	"github.com/spf13/cobra"
)

// viewCmd represents
var viewCmd = &cobra.Command{
	Use:   "view",
	Short: "read and output binary format to plain text",
	Long: `read and output binary format to plain text

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

		outFile := getFlagString(cmd, "out-file")
		showCode := getFlagBool(cmd, "show-code")
		outFasta := getFlagBool(cmd, "fasta")
		outFastq := getFlagBool(cmd, "fastq")

		outfh, gw, w, err := outStream(outFile, strings.HasSuffix(strings.ToLower(outFile), ".gz"), opt.CompressionLevel)
		checkError(err)
		defer func() {
			outfh.Flush()
			if gw != nil {
				gw.Close()
			}
			w.Close()
		}()

		var infh *bufio.Reader
		var r *os.File
		var reader *unikmer.Reader
		var kcode unikmer.KmerCode

		var quality string
		for _, file := range files {
			func() {
				infh, r, _, err = inStream(file)
				checkError(err)
				defer r.Close()

				reader, err = unikmer.NewReader(infh)
				checkError(err)

				if outFastq {
					quality = strings.Repeat("g", reader.K)
				}

				for {
					kcode, err = reader.Read()
					if err != nil {
						if err == io.EOF {
							break
						}
						checkError(err)
					}

					// outfh.WriteString(fmt.Sprintf("%s\n", kcode.Bytes())) // slower
					if outFasta {
						outfh.WriteString(fmt.Sprintf(">%d\n%s\n", kcode.Code, kcode.String()))
					} else if outFastq {
						outfh.WriteString(fmt.Sprintf(">%d\n%s\n+\n%s\n", kcode.Code, kcode.String(), quality))
					} else if showCode {
						outfh.WriteString(fmt.Sprintf("%s\t%d\n", kcode.String(), kcode.Code))
					} else {
						outfh.WriteString(kcode.String() + "\n")
					}
				}

			}()
		}
	},
}

func init() {
	RootCmd.AddCommand(viewCmd)

	viewCmd.Flags().StringP("out-file", "o", "-", `out file ("-" for stdout, suffix .gz for gzipped out)`)
	viewCmd.Flags().BoolP("show-code", "n", false, `show encoded integer along with k-mer`)
	viewCmd.Flags().BoolP("fasta", "a", false, `output in FASTA format, with encoded integer as FASTA header`)
	viewCmd.Flags().BoolP("fastq", "q", false, `output in FASTQ format, with encoded integer as FASTA header`)
}
