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
	"path/filepath"
	"runtime"
	"strings"

	"github.com/shenwei356/unikmer"
	"github.com/spf13/cobra"
)

// numCmd represents
var numCmd = &cobra.Command{
	Use:   "num",
	Short: "Quickly inspect number of k-mers in binary files",
	Long: `Quickly inspect number of k-mers in binary files

Attention:
  1. This command is designed to quickly inspect the number of k-mers in binary file,
  2. For non-sorted file, it returns '-1' unless switching on flag '-f/--force'.

`,
	Run: func(cmd *cobra.Command, args []string) {
		opt := getOptions(cmd)
		runtime.GOMAXPROCS(opt.NumCPUs)

		var err error

		if opt.Verbose {
			log.Info("checking input files ...")
		}
		files := getFileListFromArgsAndFile(cmd, args, true, "infile-list", true)
		if opt.Verbose {
			if len(files) == 1 && isStdin(files[0]) {
				log.Info("no files given, reading from stdin")
			} else {
				log.Infof("%d input file(s) given", len(files))
			}
		}

		checkFileSuffix(opt, extDataFile, files...)

		outFile := getFlagString(cmd, "out-file")
		showFile := getFlagBool(cmd, "file-name")
		basename := getFlagBool(cmd, "basename")
		force := getFlagBool(cmd, "force")

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

		for _, file := range files {
			func() {
				infh, r, _, err = inStream(file)
				checkError(err)
				defer r.Close()

				reader, err = unikmer.NewReader(infh)
				checkError(err)

				if reader.Number < 0 && force {
					var n int64
					for {
						_, _, err = reader.ReadCodeWithTaxid()
						if err != nil {
							if err == io.EOF {
								break
							}
							checkError(err)
						}

						n++
					}
					reader.Number = n
				}

				if showFile {
					if basename {
						outfh.WriteString(fmt.Sprintf("%d\t%s\n", reader.Number, filepath.Base(file)))
					} else {
						outfh.WriteString(fmt.Sprintf("%d\t%s\n", reader.Number, file))
					}
				} else {
					outfh.WriteString(fmt.Sprintf("%d\n", reader.Number))
				}
				outfh.Flush()
			}()
		}
	},
}

func init() {
	RootCmd.AddCommand(numCmd)

	numCmd.Flags().StringP("out-file", "o", "-", `out file ("-" for stdout, suffix .gz for gzipped out)`)
	numCmd.Flags().BoolP("file-name", "n", false, `show file name`)
	numCmd.Flags().BoolP("basename", "b", false, "only output basename of files")
	numCmd.Flags().BoolP("force", "f", false, "read whole file and count k-mers")
}
