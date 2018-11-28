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
	"runtime"

	"github.com/shenwei356/breader"
	"github.com/shenwei356/unikmer"
	"github.com/spf13/cobra"
)

// dumpCmd represents
var dumpCmd = &cobra.Command{
	Use:   "dump",
	Short: "convert plain k-mer text to binary format",
	Long: `convert plain k-mer text to binary format

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

		outFile := getFlagString(cmd, "out-prefix")
		unique := getFlagBool(cmd, "unique")

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

		var m map[uint64]struct{}
		if unique {
			m = make(map[uint64]struct{}, mapInitSize)
		}

		var k int = -1
		var l int
		var reader *breader.BufferedReader
		var chunk breader.Chunk
		var data interface{}
		var line string
		var kcode unikmer.KmerCode
		var ok bool
		var n int64

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

					if writer == nil {
						var mode uint32
						if opt.Compact {
							mode |= unikmer.UNIK_COMPACT
						}
						writer, err = unikmer.NewWriter(outfh, l, mode)
						checkError(err)
						defer checkError(writer.Flush())
					}

					kcode, err = unikmer.NewKmerCode([]byte(line))
					if err != nil {
						checkError(fmt.Errorf("fail to '%s': %s", line, err))
					}

					if unique {
						if _, ok = m[kcode.Code]; !ok {
							m[kcode.Code] = struct{}{}
							checkError(writer.Write(kcode))
							n++
						}
					} else {
						checkError(writer.Write(kcode))
						n++
					}
				}
			}
		}

		checkError(writer.Flush())
		if opt.Verbose {
			log.Infof("%d unique k-mers found", n)
		}
	},
}

func init() {
	RootCmd.AddCommand(dumpCmd)

	dumpCmd.Flags().StringP("out-prefix", "o", "-", `out file prefix ("-" for stdout)`)
	dumpCmd.Flags().BoolP("unique", "u", false, `remove duplicated k-mers`)
}
