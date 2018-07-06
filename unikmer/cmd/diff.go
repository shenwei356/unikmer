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

// diffCmd represents
var diffCmd = &cobra.Command{
	Use:   "diff",
	Short: "set difference of multiple binary files",
	Long: `set difference of multiple binary files

`,
	Run: func(cmd *cobra.Command, args []string) {
		opt := getOptions(cmd)
		runtime.GOMAXPROCS(opt.NumCPUs)
		files := getFileList(args)

		if len(files) < 2 {
			checkError(fmt.Errorf("at least two file should be given"))
		}

		outFile := getFlagString(cmd, "out-prefix")

		var err error

		m := make(map[uint64]bool, mapInitSize)

		var infh *xopen.Reader
		var reader *unikmer.Reader
		var kcode unikmer.KmerCode
		var k int = -1
		var firstFile = true
		var hasDiff = true
		var code uint64
		var ok bool
		for _, file := range files {
			if !isStdin(file) && !strings.HasSuffix(file, extDataFile) {
				log.Errorf("input should be stdin or %s file", extDataFile)
				return
			}

			if opt.Verbose {
				log.Infof("read kmers from %s", file)
			}

			infh, err = xopen.Ropen(file)
			checkError(err)
			defer infh.Close()

			reader, err = unikmer.NewReader(infh)
			checkError(err)

			if k == -1 {
				k = reader.K
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

				if firstFile {
					m[kcode.Code] = false
					continue
				}

				if _, ok = m[kcode.Code]; ok {
					m[kcode.Code] = true
				}
			}

			if firstFile {
				firstFile = false
				continue
			}

			// remove seen kmers
			for code = range m {
				if !m[code] {
					m[code] = false
				} else {
					delete(m, code)
				}
			}

			if len(m) == 0 {
				hasDiff = false
				break
			}
		}

		if !hasDiff {
			if opt.Verbose {
				log.Infof("no set difference found")
			}
			return
		}

		// output

		if opt.Verbose {
			log.Infof("export kmers")
		}

		if !isStdout(outFile) {
			outFile += extDataFile
		}
		outfh, err := xopen.WopenGzip(outFile)
		checkError(err)
		defer outfh.Close()

		writer := unikmer.NewWriter(outfh, k)

		for code = range m {
			writer.Write(unikmer.KmerCode{Code: code, K: k})
		}
		if opt.Verbose {
			log.Infof("set difference of %d kmers found", len(m))
		}
	},
}

func init() {
	RootCmd.AddCommand(diffCmd)

	diffCmd.Flags().StringP("out-prefix", "o", "-", `out file prefix ("-" for stdout)`)
}
