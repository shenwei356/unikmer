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

		m := make(map[uint64]struct{}, mapInitSize)

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

		k := reader.K

		// check pattern in advance
		if patternFile == "" {
			for _, query := range pattern {
				if len(query) != k {
					log.Warningf("length of query sequence (%d) != k size (%d): %s", len(query), k, query)
					return
				}
			}
		}

		for {
			kcode, err = reader.Read()
			if err != nil {
				if err == io.EOF {
					break
				}
				checkError(err)
			}

			m[kcode.Code] = struct{}{}
		}

		if opt.Verbose {
			log.Infof("finish reading kmers from %s", file)
		}

		var queries [][]byte
		var q []byte
		var ok, hit bool

		if patternFile != "" {
			var brdr *breader.BufferedReader
			brdr, err = breader.NewDefaultBufferedReader(patternFile)
			checkError(err)
			var data interface{}
			var query string
			for chunk := range brdr.Ch {
				checkError(chunk.Err)
				for _, data = range chunk.Data {
					query = data.(string)
					if query == "" {
						continue
					}
					if len(query) != k {
						log.Warningf("length of query sequence (%d) != k size (%d): %s", len(query), k, query)
						continue
					}

					query = strings.ToUpper(query)
					if degenerate {
						queries, err = extendDegenerateSeq([]byte(query))
						if err != nil {
							checkError(fmt.Errorf("extend degenerate sequence '%s': %s", query, err))
						}
					} else {
						queries = [][]byte{[]byte(query)}
					}
					for _, q = range queries {
						kcode, err = unikmer.NewKmerCode(q)
						if err != nil {
							checkError(fmt.Errorf("encoding query '%s': %s", mer, err))
						}

						_, ok = m[kcode.Code]

						if !invertMatch { //
							hit = ok
						} else {
							hit = !ok
						}

						if all {
							if hit {
								outfh.WriteString(query + "\t" + string(q) + "\n")
							}
						} else {
							if hit {
								outfh.WriteString(query + "\n")
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
					log.Warningf("length of query sequence (%d) != k size (%d): %s", len(query), k, query)
					continue
				}

				query = strings.ToUpper(query)
				if degenerate {
					queries, err = extendDegenerateSeq([]byte(query))
					if err != nil {
						checkError(fmt.Errorf("extend degenerate sequence '%s': %s", query, err))
					}
				} else {
					queries = [][]byte{[]byte(query)}
				}
				for _, q = range queries {
					kcode, err = unikmer.NewKmerCode(q)
					if err != nil {
						checkError(fmt.Errorf("encoding query '%s': %s", mer, err))
					}

					_, ok = m[kcode.Code]

					if !invertMatch { //
						hit = ok
					} else {
						hit = !ok
					}

					if all {
						if hit {
							outfh.WriteString(query + "\t" + string(q) + "\n")
						}
					} else {
						if hit {
							outfh.WriteString(query + "\n")
						}
					}
				}
			}
		}

	},
}

func init() {
	RootCmd.AddCommand(grepCmd)

	grepCmd.Flags().StringP("out-file", "o", "-", `out file ("-" for stdout, suffix .gz for gzipped out)`)

	grepCmd.Flags().StringSliceP("pattern", "p", []string{""}, `search pattern (multiple values supported. Attention: use double quotation marks for patterns containing comma, e.g., -p '"A{2,}"'))`)
	grepCmd.Flags().StringP("pattern-file", "f", "", "pattern file (one record per line)")
	grepCmd.Flags().BoolP("degenerate", "d", false, "pattern/motif contains degenerate base")
	grepCmd.Flags().BoolP("invert-match", "v", false, "invert the sense of matching, to select non-matching records")

	grepCmd.Flags().BoolP("all", "a", false, "show more information")
}
