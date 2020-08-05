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
	"sort"
	"strings"

	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/bio/seqio/fastx"
	"github.com/shenwei356/unikmer"
	"github.com/spf13/cobra"
)

// searchCmd represents
var searchCmd = &cobra.Command{
	Use:   "search",
	Short: "search sequence from index database",
	Long: `search sequence from index database

Attentions:
  0. Input format should be (gzipped) FASTA or FASTQ from file or stdin.
  1. Increase value of -j/--threads for acceleratation.
  2. Switch on -m/--use-mmap to load index files into memory to 
     accelerate searching, memory usage is roughly equal to size of index files.

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
		queryCov := getFlagFloat64(cmd, "query-cov")
		targetCov := getFlagFloat64(cmd, "target-cov")
		useMmap := getFlagBool(cmd, "use-mmap")
		nameMappingFile := getFlagString(cmd, "name-map")

		if queryCov < 0 {
			checkError(fmt.Errorf("value of -t/--query-cov should be positive"))
		}
		if targetCov < 0 {
			checkError(fmt.Errorf("value of -T/-target-cov should not be negative"))
		}

		var namesMap map[string]string
		mappingNames := nameMappingFile != ""
		if mappingNames {
			namesMap, err = readKVs(nameMappingFile, false)
			checkError(err)
			if opt.Verbose {
				log.Infof("%d pairs of name mapping values loaded", len(namesMap))
			}
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
			if useMmap {
				log.Info("loading database with mmap enabled ...")
			} else {
				log.Info("loading database ...")
			}
		}
		db, err := NewUnikIndexDB(dbDir, useMmap)
		checkError(err)
		defer func() {
			checkError(db.Close())
		}()

		if queryCov <= db.Info.FPR {
			checkError(fmt.Errorf("query coverage threshold (%f) should not small than fpr of single bloom filter of index database (%f)", queryCov, db.Info.FPR))
		}
		if opt.Verbose {
			log.Infof("db loaded: %s", db)
			log.Infof("query coverage threshold: %f", queryCov)
			log.Infof("target coverage threshold: %f", targetCov)
		}
		if mappingNames {
			var ok bool
			var _n int
			for _, name := range db.Info.Names {
				if _, ok = namesMap[name]; !ok {
					_n++
				}
			}
			if _n > 0 {
				log.Warningf("%d names are not defined in name mapping file: %s", _n, nameMappingFile)
			}
		}

		k := db.Header.K
		canonical := db.Header.Canonical

		// if !isStdout(outFile) {
		// 	outFile += ".txt"
		// }
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

					matched := db.Search(kmerList, opt.NumCPUs, queryCov, targetCov)

					targets := make([]string, 0, len(matched))
					for m := range matched {
						targets = append(targets, m)
					}
					sort.Slice(targets, func(i, j int) bool { return matched[targets[i]][0] > matched[targets[j]][0] })

					tmp := make([]string, len(targets))
					var t string
					for i, k := range targets {
						if mappingNames {
							t = namesMap[k]
						} else {
							t = k
						}
						tmp[i] = fmt.Sprintf("%s(#:%.0f, qcov:%0.4f, tcov:%0.4f)", t, matched[k][0], matched[k][1], matched[k][2])
					}

					outfh.WriteString(fmt.Sprintf("%s\t%d\t%f\t%s\n",
						record.ID, l, maxFPR(db.Info.FPR, queryCov, l), strings.Join(tmp, ", ")))
				}
			}
		}
		if opt.Verbose {
			log.Infof("done searching")
		}

	},
}

func init() {
	dbCmd.AddCommand(searchCmd)

	searchCmd.Flags().StringP("out-prefix", "o", "-", `out file prefix ("-" for stdout)`)
	searchCmd.Flags().StringP("db-dir", "d", "", `database directory created by "unikmer db index"`)
	searchCmd.Flags().Float64P("query-cov", "t", 0.6, `query coverage threshold, i.e., ratio of matched k-mers in unique k-mers of a query`)
	searchCmd.Flags().Float64P("target-cov", "T", 0, `target coverage threshold, i.e., ratio of matched k-mers in unique k-mers of a target`)
	searchCmd.Flags().BoolP("use-mmap", "m", false, `load index files into memory to accelerate searching (recommended)`)
	searchCmd.Flags().StringP("name-map", "M", "", `tabular two-column file mapping names to user-defined values`)
}
