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

	"github.com/shenwei356/breader"
	"github.com/shenwei356/unikmer"
	"github.com/shenwei356/util/pathutil"
	"github.com/spf13/cobra"
)

// grepCmd represents
var grepCmd = &cobra.Command{
	Use:   "grep",
	Short: "search Kmer from binary files",
	Long: `search Kmer from binary files

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

		if len(files) > 1 {
			checkError(fmt.Errorf("no more than one file should be given"))
		}

		checkFiles(extDataFile, files...)

		outFile := getFlagString(cmd, "out-file")
		pattern := getFlagStringSlice(cmd, "query")
		patternFile := getFlagString(cmd, "query-file")
		invertMatch := getFlagBool(cmd, "invert-match")
		degenerate := getFlagBool(cmd, "degenerate")
		all := getFlagBool(cmd, "all")

		if len(pattern) == 0 && patternFile == "" {
			checkError(fmt.Errorf("one of flags -q (--query) and -f (--query-file) needed"))
		}

		if patternFile != "" {
			var ok bool
			ok, err = pathutil.Exists(patternFile)
			if err != nil {
				checkError(fmt.Errorf("read query file: %s", err))
			}
			if !ok {
				checkError(fmt.Errorf("query file not found: %s", patternFile))
			}
		}

		file := files[0]

		m := make(map[uint64]struct{}, mapInitSize)

		if opt.Verbose {
			log.Infof("read Kmers from %s", file)
		}

		var infh *bufio.Reader
		var r *os.File
		var reader *unikmer.Reader
		var kcode unikmer.KmerCode
		var mer []byte

		infh, r, _, err = inStream(file)
		checkError(err)
		defer r.Close()

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
			log.Infof("finish reading Kmers from %s", file)
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

						if !invertMatch {
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

					if !invertMatch {
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

	grepCmd.Flags().StringSliceP("query", "q", []string{""}, `query Kmer(s) (multiple values delimted by comma supported)`)
	grepCmd.Flags().StringP("query-file", "f", "", "query file (one Kmer per line)")
	grepCmd.Flags().BoolP("degenerate", "d", false, "query Kmer contains degenerate base")
	grepCmd.Flags().BoolP("invert-match", "v", false, "invert the sense of matching, to select non-matching records")

	grepCmd.Flags().BoolP("all", "a", false, "show more information: extra column of matched Kmers")
}
