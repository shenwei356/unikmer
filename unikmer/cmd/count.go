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
	"github.com/spf13/cobra"
)

// countCmd represents
var countCmd = &cobra.Command{
	Use:   "count",
	Short: "count Kmer from FASTA/Q sequences",
	Long: `count Kmer from FASTA/Q sequences

`,
	Run: func(cmd *cobra.Command, args []string) {
		opt := getOptions(cmd)
		runtime.GOMAXPROCS(opt.NumCPUs)
		seq.ValidateSeq = false
		files := getFileList(args)

		outFile := getFlagString(cmd, "out-prefix")
		circular := getFlagBool(cmd, "circular")
		k := getFlagPositiveInt(cmd, "kmer-len")
		if k > 32 {
			checkError(fmt.Errorf("k > 32 not supported"))
		}

		canonical := getFlagBool(cmd, "canonical")

		if !isStdout(outFile) {
			outFile += extDataFile
		}
		outfh, gw, w, err := outStream(outFile, opt.Compress)
		checkError(err)
		defer func() {
			outfh.Flush()
			if gw != nil {
				gw.Close()
			}
			w.Close()
		}()

		writer := unikmer.NewWriter(outfh, k)
		writer.Compact = opt.Compact

		m := make(map[uint64]struct{}, mapInitSize)

		var sequence, kmer, preKmer []byte
		var originalLen, l, end, e int
		var record *fastx.Record
		var fastxReader *fastx.Reader
		var kcode, preKcode unikmer.KmerCode
		var first bool
		var i, j, iters int
		var ok bool
		var n int64
		for _, file := range files {
			if opt.Verbose {
				log.Infof("read sequence file: %s", file)
			}
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

				if canonical {
					iters = 1
				} else {
					iters = 2
				}

				for j = 0; j < iters; j++ {
					if j == 0 { // sequence
						sequence = record.Seq.Seq

						if opt.Verbose {
							log.Infof("process sequence: %s", record.ID)
						}
					} else { // reverse complement sequence
						sequence = record.Seq.RevComInplace().Seq

						if opt.Verbose {
							log.Infof("process reverse complement sequence: %s", record.ID)
						}
					}

					originalLen = len(record.Seq.Seq)
					l = len(sequence)

					end = l - 1
					if end < 0 {
						end = 0
					}
					first = true
					for i = 0; i <= end; i++ {
						e = i + k
						if e > originalLen {
							if circular {
								e = e - originalLen
								kmer = sequence[i:]
								kmer = append(kmer, sequence[0:e]...)
							} else {
								break
							}
						} else {
							kmer = sequence[i : i+k]
						}

						if first {
							kcode, err = unikmer.NewKmerCode(kmer)
							first = false
						} else {
							kcode, err = unikmer.NewKmerCodeMustFromFormerOne(kmer, preKmer, preKcode)
						}
						if err != nil {
							checkError(fmt.Errorf("encoding '%s': %s", kmer, err))
						}
						preKmer, preKcode = kmer, kcode

						if canonical {
							kcode = kcode.Canonical()
						}

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
	RootCmd.AddCommand(countCmd)

	countCmd.Flags().StringP("out-prefix", "o", "-", `out file prefix ("-" for stdout)`)
	countCmd.Flags().IntP("kmer-len", "k", 0, "Kmer length")
	countCmd.Flags().BoolP("circular", "", false, "circular genome")
	countCmd.Flags().BoolP("canonical", "", false, "only keep the canonical Kmers")
}
