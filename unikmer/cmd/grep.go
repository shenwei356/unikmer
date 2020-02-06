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
	"bufio"
	"fmt"
	"io"
	"os"
	"runtime"
	"strings"

	"github.com/shenwei356/breader"
	"github.com/shenwei356/unikmer"
	"github.com/spf13/cobra"
)

// grepCmd represents
var grepCmd = &cobra.Command{
	Use:   "grep",
	Short: "search k-mers from binary files",
	Long: `search k-mers from binary files

`,
	Run: func(cmd *cobra.Command, args []string) {
		opt := getOptions(cmd)
		runtime.GOMAXPROCS(opt.NumCPUs)

		var err error

		files := getFileListFromArgsAndFile(cmd, args, true, "infile-list", true)

		if len(files) > 1 {
			checkError(fmt.Errorf("no more than one file should be given"))
		}

		checkFileSuffix(extDataFile, files...)

		outFile := getFlagString(cmd, "out-file")
		pattern := getFlagStringSlice(cmd, "query")
		queryFile := getFlagString(cmd, "query-file")
		invertMatch := getFlagBool(cmd, "invert-match")
		degenerate := getFlagBool(cmd, "degenerate")
		canonicalSensitive := getFlagBool(cmd, "canonical")
		// all := getFlagBool(cmd, "all")

		if len(pattern) == 0 && queryFile == "" {
			checkError(fmt.Errorf("one of flags -q (--query) and -f (--query-file) needed"))
		}

		m := make(map[uint64]struct{}, mapInitSize)

		k := -1
		queryList := make([]string, 0, 1000)
		for _, query := range pattern {
			if query == "" {
				continue
			}
			if k == -1 {
				k = len(query)
			}
			if len(query) != k {
				checkError(fmt.Errorf("length of query sequence are inconsistent: (%d) != (%d): %s", len(query), k, query))
			}
			queryList = append(queryList, query)
		}
		if queryFile != "" {
			if opt.Verbose {
				log.Infof("load queries from file: %s", queryFile)
			}
			var brdr *breader.BufferedReader
			brdr, err = breader.NewDefaultBufferedReader(queryFile)
			checkError(err)
			var data interface{}
			var query string
			for chunk := range brdr.Ch {
				checkError(chunk.Err)
				for _, data = range chunk.Data {
					query = data.(string)
					if k == -1 {
						k = len(query)
					}
					if len(query) != k {
						checkError(fmt.Errorf("length of query sequence are inconsistent: (%d) != (%d): %s", len(query), k, query))

					}
					queryList = append(queryList, query)
				}
			}
		}

		var kcode unikmer.KmerCode
		var mer []byte

		var queries [][]byte
		var q []byte
		for _, query := range queryList {
			if degenerate {
				queries, err = extendDegenerateSeq([]byte(query))
				if err != nil {
					checkError(fmt.Errorf("fail to extend degenerate sequence '%s': %s", query, err))
				}
			} else {
				queries = [][]byte{[]byte(query)}
			}

			for _, q = range queries {
				kcode, err = unikmer.NewKmerCode(q)
				if err != nil {
					checkError(fmt.Errorf("fail to encode query '%s': %s", mer, err))
				}
				m[kcode.Code] = struct{}{}
			}
		}

		if opt.Verbose {
			log.Infof("%d queries loaded", len(m))
			if canonicalSensitive {
				log.Infof("will check reverse complement sequence of canonical k-mers")
			}
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

		var n int64
		var infh *bufio.Reader
		var r *os.File
		var reader *unikmer.Reader
		var canonical bool
		var firstFile = true
		var flag int
		var nfiles = len(files)
		var ok, hit bool
		for i, file := range files {
			if !firstFile && file == files[0] {
				continue
			}

			if opt.Verbose {
				log.Infof("processing file (%d/%d): %s", i+1, nfiles, file)
			}

			flag = func() int {
				infh, r, _, err = inStream(file)
				checkError(err)
				defer r.Close()

				reader, err = unikmer.NewReader(infh)
				checkError(err)

				if k != reader.K {
					checkError(fmt.Errorf("K (%d) of binary file '%s' not equal to previous K (%d)", reader.K, file, k))
				}
				canonical = reader.Flag&unikmer.UNIK_CANONICAL > 0
				if opt.Verbose {
					log.Infof("file (%d/%d): %s: canonical: %v", i+1, nfiles, file, canonical)
				}

				for {
					kcode, err = reader.Read()
					if err != nil {
						if err == io.EOF {
							break
						}
						checkError(err)
					}

					_, ok = m[kcode.Code]

					if canonicalSensitive && canonical && !ok {
						_, ok = m[kcode.RevComp().Code]
					}

					if !invertMatch {
						hit = ok
					} else {
						hit = !ok
					}

					if hit {
						outfh.WriteString(kcode.String() + "\n")
						outfh.Flush()
						n++
					}
				}

				return flagContinue
			}()

			if flag == flagReturn {
				return
			} else if flag == flagBreak {
				break
			}
		}

		if opt.Verbose {
			log.Infof("%d k-mer(s) found", n)
		}

	},
}

func init() {
	RootCmd.AddCommand(grepCmd)

	grepCmd.Flags().StringP("out-file", "o", "-", `out file ("-" for stdout, suffix .gz for gzipped out)`)

	grepCmd.Flags().StringSliceP("query", "q", []string{""}, `query k-mers (multiple values delimted by comma supported)`)
	grepCmd.Flags().StringP("query-file", "f", "", "query file (one k-mer per line)")
	grepCmd.Flags().BoolP("degenerate", "d", false, "query k-mers contains degenerate base")
	grepCmd.Flags().BoolP("invert-match", "v", false, "invert the sense of matching, to select non-matching records")
	grepCmd.Flags().BoolP("canonical", "K", false, "check reverse complement sequence of canonical k-mers")
	// grepCmd.Flags().BoolP("all", "a", false, "show more information: extra column of matched k-mers")
}
