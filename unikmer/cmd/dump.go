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

		outFile := getFlagString(cmd, "out-prefix")
		unique := getFlagBool(cmd, "unique")
		canonical := getFlagBool(cmd, "canonical")
		canonicalOnly := getFlagBool(cmd, "canonical-only")
		sortedKmers := getFlagBool(cmd, "sorted")

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
		var kcode, kcodeC unikmer.KmerCode
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
						if sortedKmers {
							mode |= unikmer.UNIK_SORTED
						} else if opt.Compact {
							mode |= unikmer.UNIK_COMPACT
						}
						if canonical || canonicalOnly {
							mode |= unikmer.UNIK_CANONICAL
						}
						writer, err = unikmer.NewWriter(outfh, l, mode)
						checkError(err)
					}

					kcode, err = unikmer.NewKmerCode([]byte(line))
					if err != nil {
						checkError(fmt.Errorf("fail to encode '%s': %s", line, err))
					}

					if canonicalOnly {
						kcodeC = kcode.Canonical()
						if kcode.Code != kcodeC.Code {
							continue
						}
						kcode = kcodeC
					} else if canonical {
						kcode = kcode.Canonical()
					}

					if unique {
						if _, ok = m[kcode.Code]; !ok {
							m[kcode.Code] = struct{}{}
							writer.WriteCode(kcode.Code)
							n++
						}
					} else {
						writer.WriteCode(kcode.Code)
						n++
					}
				}
			}
		}

		checkError(writer.Flush())
		if opt.Verbose {
			log.Infof("%d unique k-mers saved to %s", n, outFile)
		}
	},
}

func init() {
	RootCmd.AddCommand(dumpCmd)

	dumpCmd.Flags().StringP("out-prefix", "o", "-", `out file prefix ("-" for stdout)`)
	dumpCmd.Flags().BoolP("unique", "u", false, `remove duplicated k-mers`)
	dumpCmd.Flags().BoolP("canonical", "K", false, "save the canonical k-mers")
	dumpCmd.Flags().BoolP("canonical-only", "k", false, "only save the canonical k-mers. This option overides -K/--canonical")
	dumpCmd.Flags().BoolP("sorted", "s", false, "input k-mers are sorted")
}
