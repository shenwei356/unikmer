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
	"github.com/shenwei356/xopen"
	"github.com/spf13/cobra"
)

// dumpCmd represents
var dumpCmd = &cobra.Command{
	Use:   "dump",
	Short: "convert plain Kmer text to binary format",
	Long: `convert plain Kmer text to binary format

`,
	Run: func(cmd *cobra.Command, args []string) {
		opt := getOptions(cmd)
		runtime.GOMAXPROCS(opt.NumCPUs)
		files := getFileList(args)

		outFile := getFlagString(cmd, "out-prefix")
		noDedup := getFlagBool(cmd, "no-dedup")

		if !isStdout(outFile) {
			outFile += extDataFile
		}
		outfh, err := xopen.WopenGzip(outFile)
		checkError(err)
		defer outfh.Close()

		var writer *unikmer.Writer

		var m map[uint64]struct{}
		if !noDedup {
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
						checkError(fmt.Errorf("Kmer length mismatch, previous: %d, current: %d. %s", k, l, line))
					}

					if writer == nil {
						writer = unikmer.NewWriter(outfh, l)
					}

					kcode, err = unikmer.NewKmerCode([]byte(line))
					if err != nil {
						checkError(fmt.Errorf("encoding '%s': %s", line, err))
					}

					if noDedup {
						checkError(writer.Write(kcode))
						n++
					} else {
						if _, ok = m[kcode.Code]; !ok {
							m[kcode.Code] = struct{}{}
							checkError(writer.Write(kcode))
							n++
						}
					}
				}
			}
		}

		if opt.Verbose {
			log.Infof("%d unique Kmers found", n)
		}
	},
}

func init() {
	RootCmd.AddCommand(dumpCmd)

	dumpCmd.Flags().StringP("out-prefix", "o", "-", `out file prefix ("-" for stdout)`)
	dumpCmd.Flags().BoolP("no-dedup", "D", false, `do not deduplicate kmers, this can save some time and memory`)
}
