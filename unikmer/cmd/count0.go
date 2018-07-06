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

	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/bio/seqio/fastx"
	"github.com/shenwei356/unikmer"
	"github.com/shenwei356/xopen"
	"github.com/spf13/cobra"
)

// mcountCmd represents
var mcountCmd = &cobra.Command{
	Use:   "mcount",
	Short: "count kmer from FASTA/Q sequences",
	Long: `count kmer from FASTA/Q sequences

`,
	Run: func(cmd *cobra.Command, args []string) {
		opt := getOptions(cmd)
		runtime.GOMAXPROCS(opt.NumCPUs)
		seq.ValidateSeq = false
		files := getFileList(args)

		outFile := getFlagString(cmd, "out-prefix")
		circular := getFlagBool(cmd, "circular-genome")
		k := getFlagPositiveInt(cmd, "kmer-len")
		if k > 32 {
			checkError(fmt.Errorf("k > 32 not supported"))
		}
		// hint := getFlagPositiveInt(cmd, "esti-kmer-num")

		if !isStdout(outFile) {
			outFile += extDataFile
		}
		outfh, err := xopen.WopenGzip(outFile)
		checkError(err)
		defer outfh.Close()

		writer := unikmer.NewWriter(outfh, k)

		// estimate according to value of  k
		// sbf := boom.NewScalableBloomFilter(uint(hint), 0.01, 0.8)
		m := make(map[uint64]struct{}, 100000)

		var sequence, mer []byte // , merRC []byte
		var originalLen, l, end, e int
		var record *fastx.Record
		var fastxReader *fastx.Reader
		var kcode unikmer.KmerCode
		var i int
		var ok bool
		for _, file := range files {
			fastxReader, err = fastx.NewDefaultReader(file)
			checkError(err)
			for {
				record, err = fastxReader.Read()
				if err != nil {
					if err == io.EOF {
						break
					}
					checkError(err)
					break
				}

				originalLen = len(record.Seq.Seq)
				sequence = record.Seq.Seq
				l = len(sequence)
				end = l - 1
				if end < 0 {
					end = 0
				}
				for i = 0; i <= end; i++ {
					e = i + k
					if e > originalLen {
						if circular {
							e = e - originalLen
							mer = sequence[i:]
							mer = append(mer, sequence[0:e]...)
						} else {
							break
						}
					} else {
						mer = sequence[i : i+k]
					}

					kcode, err = unikmer.NewKmerCode(mer)
					if err != nil {
						checkError(fmt.Errorf("encoding '%s': %s", mer, err))
					}

					// mer = kcode.Bytes()
					if _, ok = m[kcode.Code]; !ok {
						// if !sbf.Test(mer) {
						// sbf.Add(mer)
						m[kcode.Code] = struct{}{}
						checkError(writer.Write(kcode))
					}
				}

				// ------------------

				sequence = record.Seq.RevComInplace().Seq
				l = len(sequence)
				end = l - 1
				if end < 0 {
					end = 0
				}
				for i = 0; i <= end; i++ {
					e = i + k
					if e > originalLen {
						if circular {
							e = e - originalLen
							mer = sequence[i:]
							mer = append(mer, sequence[0:e]...)
						} else {
							break
						}
					} else {
						mer = sequence[i : i+k]
					}

					kcode, err = unikmer.NewKmerCode(mer)
					if err != nil {
						checkError(fmt.Errorf("encoding '%s': %s", mer, err))
					}

					// mer = kcode.Bytes()
					if _, ok = m[kcode.Code]; !ok {
						// if !sbf.Test(mer) {
						// sbf.Add(mer)
						m[kcode.Code] = struct{}{}
						checkError(writer.Write(kcode))
					}
				}

			}
		}

	},
}

func init() {
	RootCmd.AddCommand(mcountCmd)

	mcountCmd.Flags().StringP("out-prefix", "o", "-", `out file prefix ("-" for stdout)`)
	mcountCmd.Flags().IntP("kmer-len", "k", 0, "kmer length")
	// mcountCmd.Flags().IntP("esti-kmer-num", "n", 100000000, "estimated kmer num length (for initializing Bloom Filter)")
	mcountCmd.Flags().BoolP("circular-genome", "c", false, "circular genome")
}
