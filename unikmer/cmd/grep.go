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
	"compress/gzip"
	"fmt"
	"io"
	"os"
	"runtime"
	"sync"

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
		// invertMatch := getFlagBool(cmd, "invert-match")
		degenerate := getFlagBool(cmd, "degenerate")

		if len(pattern) == 0 && patternFile == "" {
			checkError(fmt.Errorf("one of flags -p (--pattern) and -f (--pattern-file) needed"))
		}
		ok, err := pathutil.Exists(patternFile)
		if patternFile != "" && !ok {
			checkError(fmt.Errorf("pattern file not found: %s", patternFile))
		}

		outfh, err := xopen.Wopen(outFile)
		checkError(err)
		defer outfh.Close()

		file := files[0]

		fileSBF, fileIBF := file+extSBF, file+extIBF
		fileSBFExists, err := pathutil.Exists(fileSBF)
		checkError(err)
		fileIBFExists, err := pathutil.Exists(fileIBF)
		checkError(err)

		var sbf *boom.ScalableBloomFilter
		var ibf *boom.InverseBloomFilter

		var headerSBF, headerIBF unikmer.Header

		var k int

		if !isStdin(file) && fileSBFExists && fileIBFExists { // read index file
			var wg sync.WaitGroup

			wg.Add(1)
			go func() {
				defer wg.Done()

				r, err := os.Open(fileSBF)
				checkError(err)
				defer r.Close()

				br := bufio.NewReader(r)

				infh, err := gzip.NewReader(br)
				checkError(err)
				defer infh.Close()

				headerSBF, err = readHeader(infh)
				if err != nil {
					checkError(fmt.Errorf("read .sbf file '%s': %s", fileSBF, err))
				}

				sbf = &boom.ScalableBloomFilter{}
				_, err = sbf.ReadFrom(infh)
				checkError(err)
			}()

			wg.Add(1)
			go func() {
				defer wg.Done()

				r, err := os.Open(fileIBF)
				checkError(err)
				defer r.Close()

				br := bufio.NewReader(r)

				infh, err := gzip.NewReader(br)
				checkError(err)
				defer infh.Close()

				headerIBF, err = readHeader(infh)
				if err != nil {
					checkError(fmt.Errorf("read .ibf file '%s': %s", fileIBF, err))
				}

				ibf = &boom.InverseBloomFilter{}
				_, err = ibf.ReadFrom(infh)
				checkError(err)
			}()

			wg.Wait()

			if headerIBF.Version != headerSBF.Version {
				checkError(fmt.Errorf("version mismatch: %s (.sbf) != %s (.ibf)", headerSBF.Version, headerIBF.Version))
			}
			if headerIBF.K != headerSBF.K {
				checkError(fmt.Errorf("k size mismatch: %d (.sbf) != %d (.ibf)", headerSBF.K, headerIBF.K))
			}

			k = headerIBF.K

		} else { // read binary file
			sbf = boom.NewScalableBloomFilter(uint(hint), 0.01, 0.8)
			ibf = boom.NewInverseBloomFilter(uint(hint / 5))

			var infh *xopen.Reader
			var reader *unikmer.Reader
			var kcode unikmer.KmerCode
			var mer []byte
			var err error

			infh, err = xopen.Ropen(file)
			checkError(err)
			defer infh.Close()

			reader, err = unikmer.NewReader(infh)
			checkError(err)

			sbf = boom.NewScalableBloomFilter(uint(hint), 0.01, 0.8)
			ibf = boom.NewInverseBloomFilter(uint(hint / 5))

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

			k = reader.K
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
					query = data.(string)
					if query == "" {
						continue
					}
					if len(query) != k {
						log.Errorf("length of query sequence (%d) != k size (%d)", len(query), k)
					}

					if degenerate {
						queries = make([]string, 0, 4)
					} else {
						queries = []string{query}
					}

					for _, q = range queries {
						qb = []byte(q)
						inSBF = sbf.Test(qb)
						inIBF = ibf.Test(qb)

						fmt.Printf("%s\t%v\t%v\n", q, inSBF, inIBF)

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
}
