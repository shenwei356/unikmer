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

	"github.com/shenwei356/breader"
	"github.com/shenwei356/unikmer"
	"github.com/shenwei356/util/pathutil"
	"github.com/shenwei356/xopen"
	"github.com/spf13/cobra"
	boom "github.com/tylertreat/BoomFilters"
)

// grepCmd represents
var grepCmd = &cobra.Command{
	Use:   "grep",
	Short: "search kmer from binary file",
	Long: `search kmer from binary file

`,
	Run: func(cmd *cobra.Command, args []string) {
		opt := getOptions(cmd)
		runtime.GOMAXPROCS(opt.NumCPUs)
		files := getFileList(args)

		if len(files) > 1 {
			checkError(fmt.Errorf("no more than one file should be given"))
		}

		outFile := getFlagString(cmd, "out-file")
		hint := getFlagPositiveInt(cmd, "esti-kmer-num")
		pattern := getFlagStringSlice(cmd, "pattern")
		patternFile := getFlagString(cmd, "pattern-file")
		invertMatch := getFlagBool(cmd, "invert-match")
		degenerate := getFlagBool(cmd, "degenerate")
		all := getFlagBool(cmd, "all")

		if len(pattern) == 0 && patternFile == "" {
			checkError(fmt.Errorf("one of flags -p (--pattern) and -f (--pattern-file) needed"))
		}

		var err error

		if patternFile != "" {
			var ok bool
			ok, err = pathutil.Exists(patternFile)
			if !ok {
				checkError(fmt.Errorf("read pattern file: %s", err))
			}
			if !ok {
				checkError(fmt.Errorf("pattern file not found: %s", patternFile))
			}
		}

		outfh, err := xopen.Wopen(outFile)
		checkError(err)
		defer outfh.Close()

		file := files[0]

		if !isStdin(file) && !strings.HasSuffix(file, extDataFile) {
			log.Errorf("input should be stdin or %s file", extDataFile)
			return
		}

		var sbf *boom.ScalableBloomFilter
		var ibf *boom.InverseBloomFilter

		sbf = boom.NewScalableBloomFilter(uint(hint), 0.01, 0.8)
		ibf = boom.NewInverseBloomFilter(uint(hint / 5))

		if opt.Verbose {
			log.Infof("read kmers from %s", file)
		}

		var infh *xopen.Reader
		var reader *unikmer.Reader
		var kcode unikmer.KmerCode
		var mer []byte
		infh, err = xopen.Ropen(file)
		checkError(err)
		defer infh.Close()

		reader, err = unikmer.NewReader(infh)
		checkError(err)

		for {
			kcode, err = reader.Read()
			if err != nil {
				if err == io.EOF {
					break
				}
				checkError(err)
			}

			mer = kcode.Bytes()
			sbf.Add(mer)
			ibf.Add(mer)
		}

		k := reader.K

		if opt.Verbose {
			log.Infof("finish reading kmers from %s", file)
		}

		var queries []string
		var q string
		var qb []byte
		var inSBF, inIBF bool

		if patternFile != "" {
			reader, err := breader.NewDefaultBufferedReader(patternFile)
			checkError(err)
			var data interface{}
			var query string
			for chunk := range reader.Ch {
				checkError(chunk.Err)
				for _, data = range chunk.Data {
					query = strings.ToUpper(data.(string))
					if query == "" {
						continue
					}
					if len(query) != k {
						log.Warningf("length of query sequence (%d) != k size (%d)", len(query), k)
					}

					if degenerate {
						queries = make([]string, 0, 4)
						// todo
					} else {
						queries = []string{query}
					}

					for _, q = range queries {
						qb = []byte(q)
						inSBF = sbf.Test(qb)
						inIBF = ibf.Test(qb)

						fmt.Printf("%s\t%v\t%v\n", q, inSBF, inIBF)

						// - Scalable Bloom Filter: false-negative == 0, 不在就真不在
						// - Inverse Bloom Filter : false positive == 0，在就真的在
						if !invertMatch { //
							if inIBF {
								if all {
									outfh.WriteString(query + "\t" + q + "\n")
								} else {
									outfh.WriteString(query + "\n")
								}
							}
						} else {
							if !inSBF {
								if all {
									outfh.WriteString(query + "\t" + q + "\n")
								} else {
									outfh.WriteString(query + "\n")
								}
							}
						}

					}

				}
			}
		} else {
			for _, query := range pattern {
				if query == "" {
					continue
				}
				if len(query) != k {
					log.Errorf("length of query sequence (%d) != k size (%d)", len(query), k)
				}

			}
		}

	},
}

func init() {
	RootCmd.AddCommand(grepCmd)

	grepCmd.Flags().StringP("out-file", "o", "-", `out file ("-" for stdout, suffix .gz for gzipped out)`)
	grepCmd.Flags().IntP("esti-kmer-num", "n", 100000000, "estimated kmer num length (for initializing Bloom Filter)")

	grepCmd.Flags().StringSliceP("pattern", "p", []string{""}, `search pattern (multiple values supported. Attention: use double quotation marks for patterns containing comma, e.g., -p '"A{2,}"'))`)
	grepCmd.Flags().StringP("pattern-file", "f", "", "pattern file (one record per line)")
	grepCmd.Flags().BoolP("degenerate", "d", false, "pattern/motif contains degenerate base")
	grepCmd.Flags().BoolP("invert-match", "v", false, "invert the sense of matching, to select non-matching records")

	grepCmd.Flags().BoolP("all", "a", false, "show more information")
}
