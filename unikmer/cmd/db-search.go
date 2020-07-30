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
	"io"
	"runtime"

	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/bio/seqio/fastx"
	"github.com/shenwei356/unikmer"
	"github.com/spf13/cobra"
)

// searchCmd represents
var searchCmd = &cobra.Command{
	Use:   "search",
	Short: "search sequence",
	Long: `search sequence

`,
	Run: func(cmd *cobra.Command, args []string) {
		opt := getOptions(cmd)
		runtime.GOMAXPROCS(opt.NumCPUs)
		seq.ValidateSeq = false

		var err error

		dbDir := getFlagString(cmd, "db-dir")
		if dbDir == "" {
			checkError(fmt.Errorf("flag -d/--db-dir needed"))
		}
		outFile := getFlagString(cmd, "out-prefix")

		moreVerbose := getFlagBool(cmd, "more-verbose")

		if moreVerbose {
			opt.Verbose = true
		}

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

		if opt.Verbose {
			log.Info("loading database ...")
		}
		db, err := NewUnikIndexDB(dbDir)
		checkError(err)
		defer func() {
			checkError(db.Close())
		}()

		if opt.Verbose {
			log.Infof("db loaded: %s", db)
		}
		k := db.Header.K
		canonical := db.Header.Canonical

		if !isStdout(outFile) {
			outFile += ".txt"
		}
		outfh, gw, w, err := outStream(outFile, false, opt.CompressionLevel)
		checkError(err)
		defer func() {
			outfh.Flush()
			if gw != nil {
				gw.Close()
			}
			w.Close()
		}()

		var sequence, kmer, preKmer []byte
		var originalLen, l, end, e int
		var record *fastx.Record
		var fastxReader *fastx.Reader
		var kcode, preKcode unikmer.KmerCode
		var first bool
		var i, j, iters int
		var nseq int64
		for _, file := range files {
			if opt.Verbose {
				log.Infof("reading sequence file: %s", file)
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

				nseq++
				if opt.Verbose && moreVerbose {
					log.Infof("searching sequence #%d: %s", nseq, record.ID)
				}

				if canonical {
					iters = 1
				} else {
					iters = 2
				}

				for j = 0; j < iters; j++ {
					if j == 0 { // sequence
						sequence = record.Seq.Seq
					} else { // reverse complement sequence
						sequence = record.Seq.RevComInplace().Seq
					}

					kmers := make(map[uint64]int, 128)
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
							break
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
							checkError(fmt.Errorf("fail to encode '%s': %s", kmer, err))
						}
						preKmer, preKcode = kmer, kcode

						if canonical {
							kcode = kcode.Canonical()
						}

						kmers[kcode.Code]++
					}

					kmerList := make([]uint64, 0, len(kmers))
					for code := range kmers {
						kmerList = append(kmerList, code)
					}

					result := db.Search(kmerList, opt.NumCPUs)
					outfh.WriteString(fmt.Sprintf("%s: %s", record.ID, result))
				}
			}
		}

	},
}

func init() {
	dbCmd.AddCommand(searchCmd)

	searchCmd.Flags().StringP("out-prefix", "o", "-", `out file prefix ("-" for stdout)`)
	searchCmd.Flags().StringP("db-dir", "d", "", `database directory created by "unikmer db index"`)
	searchCmd.Flags().BoolP("canonical", "K", false, "only keep the canonical k-mers")
	searchCmd.Flags().BoolP("more-verbose", "V", false, `print extra verbose information`)

}
