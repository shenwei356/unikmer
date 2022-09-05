// Copyright © 2018-2021 Wei Shen <shenwei356@gmail.com>
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
	"strings"

	"github.com/pkg/errors"
	"github.com/shenwei356/breader"
	"github.com/shenwei356/kmers"

	"github.com/spf13/cobra"
	"github.com/will-rowe/nthash"
)

var encodeCmd = &cobra.Command{
	Use:   "encode",
	Short: "Encode plain k-mer text to integer",
	Long: `Encode plain k-mer text to integer

`,
	Run: func(cmd *cobra.Command, args []string) {
		opt := getOptions(cmd)

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

		outFile := getFlagString(cmd, "out-file")
		all := getFlagBool(cmd, "all")
		canonical := getFlagBool(cmd, "canonical")
		hashed := getFlagBool(cmd, "hash")

		outfh, gw, w, err := outStream(outFile, strings.HasSuffix(strings.ToLower(outFile), ".gz"), opt.CompressionLevel)
		checkError(err)
		defer func() {
			outfh.Flush()
			if gw != nil {
				gw.Close()
			}
			w.Close()
		}()

		var k int = -1
		var l int
		var reader *breader.BufferedReader
		var chunk breader.Chunk
		var data interface{}
		var line string
		var linebytes []byte
		var kcode kmers.KmerCode
		var hasher *nthash.NTHi
		var hash uint64
		for _, file := range files {
			reader, err = breader.NewDefaultBufferedReader(file)
			checkError(errors.Wrap(err, file))

			for chunk = range reader.Ch {
				checkError(chunk.Err)
				for _, data = range chunk.Data {
					line = data.(string)
					l = len(line)

					if l == 0 {
						continue
					} else if k == -1 {
						k = l
						if k > 64 {
							checkError(fmt.Errorf("k-mer size (%d) should be <=64", k))
						}
					} else if l != k {
						checkError(fmt.Errorf("K-mer length mismatch, previous: %d, current: %d. %s", k, l, line))
					}

					if hashed {
						linebytes = []byte(line)
						hasher, err = nthash.NewHasher(&linebytes, uint(k))
						checkError(errors.Wrap(err, line))
						// for hash = range hasher.Hash(canonical) {
						hash, _ = hasher.Next(canonical)
						outfh.WriteString(fmt.Sprintf("%d\n", hash))

						continue
					}

					kcode, err = kmers.NewKmerCode([]byte(line))
					if err != nil {
						checkError(fmt.Errorf("fail to encode '%s': %s", line, err))
					}
					if canonical {
						kcode = kcode.Canonical()
					}

					if all {
						outfh.WriteString(fmt.Sprintf("%s\t%s\t%d\t%s\n", line, kcode.String(), kcode.Code, kcode.BitsString()))
					} else {
						outfh.WriteString(fmt.Sprintf("%d\n", kcode.Code))
					}
				}
			}
		}

	},
}

func init() {
	RootCmd.AddCommand(encodeCmd)

	encodeCmd.Flags().StringP("out-file", "o", "-", `out file ("-" for stdout, suffix .gz for gzipped out)`)
	encodeCmd.Flags().BoolP("all", "a", false, `output all data: orginial k-mer, parsed k-mer, encoded integer, encode bits`)
	encodeCmd.Flags().BoolP("canonical", "K", false, "keep the canonical k-mers")
	encodeCmd.Flags().BoolP("hash", "H", false, `save hash of k-mer, automatically on for k>32`)
}
