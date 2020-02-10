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
	"bufio"
	"fmt"
	"io"
	"os"
	"runtime"
	"strings"

	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/bio/seqio/fastx"
	"github.com/shenwei356/unikmer"
	"github.com/spf13/cobra"
)

// uniqsCmd represents
var uniqsCmd = &cobra.Command{
	Use:   "uniqs",
	Short: "mapping k-mers back to genome and find unique subsequences",
	Long: `mapping k-mers back to genome and find unique subsequences

Attention:
  1. default output is in BED3 format, with left-closed and right-open
     0-based interval
`,
	Run: func(cmd *cobra.Command, args []string) {
		opt := getOptions(cmd)
		runtime.GOMAXPROCS(opt.NumCPUs)
		seq.ValidateSeq = false

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

		checkFileSuffix(extDataFile, files...)

		outFile := getFlagString(cmd, "out-prefix")
		circular := getFlagBool(cmd, "circular")

		genomeFile := getFlagNonEmptyString(cmd, "genome")

		minLen := getFlagPositiveInt(cmd, "min-len")
		mMapped := getFlagBool(cmd, "allow-muliple-mapped-kmer")
		outputFASTA := getFlagBool(cmd, "output-fasta")
		maxContNonUniqKmers := getFlagNonNegativeInt(cmd, "max-cont-non-uniq-kmers")
		maxContNonUniqKmersNum := getFlagNonNegativeInt(cmd, "max-num-cont-non-uniq-kmers")

		if maxContNonUniqKmersNum > 0 && maxContNonUniqKmers == 0 {
			log.Warningf("-X/--max-num-cont-non-uniq-kmers %d is ignored becaue value of -x/--max-cont-non-uniq-kmers is 0", maxContNonUniqKmersNum)
		}
		if maxContNonUniqKmers > 0 && maxContNonUniqKmersNum == 0 {
			checkError(fmt.Errorf("value of -X/--max-num-cont-non-uniq-kmers should be > 0 when value of -x/--max-cont-non-uniq-kmers is > 0"))
		}

		m := make(map[uint64]struct{}, mapInitSize)

		// debug := bufio.NewWriterSize(os.Stdout, os.Getpagesize())
		// defer debug.Flush()

		// -----------------------------------------------------------------------

		var k int = -1
		var canonical bool

		var infh *bufio.Reader
		var r *os.File
		var reader *unikmer.Reader
		var kcode unikmer.KmerCode
		var nfiles = len(files)
		for i, file := range files {
			if opt.Verbose {
				log.Infof("reading file (%d/%d): %s", i+1, nfiles, file)
			}
			func() {
				infh, r, _, err = inStream(file)
				checkError(err)
				defer r.Close()

				reader, err = unikmer.NewReader(infh)
				checkError(err)

				if k == -1 {
					k = reader.K
					canonical = reader.Flag&unikmer.UNIK_CANONICAL > 0
					if opt.Verbose {
						if canonical {
							log.Infof("flag of canonical is on")
						} else {
							log.Infof("flag of canonical is off")
						}
					}
				} else if k != reader.K {
					checkError(fmt.Errorf("K (%d) of binary file '%s' not equal to previous K (%d)", reader.K, file, k))
				} else if (reader.Flag&unikmer.UNIK_CANONICAL > 0) != canonical {
					checkError(fmt.Errorf(`'canonical' flags not consistent, please check with "unikmer stats"`))
				}

				if canonical {
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
				} else {
					for {
						kcode, err = reader.Read()
						if err != nil {
							if err == io.EOF {
								break
							}
							checkError(err)
						}

						m[kcode.Canonical().Code] = struct{}{}
					}
				}
			}()
		}

		if opt.Verbose {
			log.Infof("%d k-mers loaded", len(m))
		}

		// -----------------------------------------------------------------------
		var m2 map[uint64]bool

		var sequence, kmer, preKmer []byte
		var originalLen, l, end, e int
		var record *fastx.Record
		var fastxReader *fastx.Reader
		var preKcode unikmer.KmerCode
		var first bool
		var i int
		var ok bool
		var multipleMapped bool

		if !mMapped {
			m2 = make(map[uint64]bool, mapInitSize)
			if opt.Verbose {
				log.Infof("pre-reading genome file: %s", genomeFile)
			}
			fastxReader, err = fastx.NewDefaultReader(genomeFile)
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

				sequence = record.Seq.Seq

				if opt.Verbose {
					log.Infof("processing sequence: %s", record.ID)
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
						checkError(fmt.Errorf("fail to encode '%s': %s", kmer, err))
					}
					preKmer, preKcode = kmer, kcode

					kcode = kcode.Canonical()

					if multipleMapped, ok = m2[kcode.Code]; !ok {
						m2[kcode.Code] = false
					} else if !multipleMapped {
						m2[kcode.Code] = true
					}
				}
			}
			if opt.Verbose {
				log.Infof("finished pre-reading genome file: %s", genomeFile)
			}

			if opt.Verbose {
				log.Infof("%d k-mers loaded from genome", len(m2))
			}
			for code, flag := range m2 {
				if !flag {
					delete(m2, code)
				}
			}
			if opt.Verbose {
				log.Infof("%d k-mers in genome are multiple mapped", len(m2))
			}
		}

		// -----------------------------------------------------------------------

		outfh, gw, w, err := outStream(outFile, strings.HasSuffix(strings.ToLower(outFile), ".gz"), opt.CompressionLevel)
		checkError(err)
		defer func() {
			outfh.Flush()
			if gw != nil {
				gw.Close()
			}
			w.Close()
		}()

		var c, start, nonUniqs, nonUniqsNum, lastNonUniqsNum, lastmatch int
		var flag bool = true
		if opt.Verbose {
			log.Infof("reading genome file: %s", genomeFile)
		}
		fastxReader, err = fastx.NewDefaultReader(genomeFile)
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

			sequence = record.Seq.Seq

			if opt.Verbose {
				log.Infof("processinig sequence: %s", record.ID)
			}

			originalLen = len(record.Seq.Seq)
			l = len(sequence)

			end = l - 1
			if end < 0 {
				end = 0
			}

			c = 0
			start = -1
			nonUniqs = 0
			nonUniqsNum = 0

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

				kcode = kcode.Canonical()

				if _, ok = m[kcode.Code]; ok {
					nonUniqs = 0
					if !mMapped {
						if multipleMapped, ok = m2[kcode.Code]; ok && multipleMapped {
							if lastNonUniqsNum <= maxContNonUniqKmersNum &&
								start >= 0 && lastmatch-start+k >= minLen {
								if outputFASTA {
									outfh.WriteString(fmt.Sprintf(">%s:%d-%d\n%s\n", record.ID, start+1, lastmatch+k,
										record.Seq.SubSeq(start+1, lastmatch+k).FormatSeq(60)))
								} else {
									outfh.WriteString(fmt.Sprintf("%s\t%d\t%d\n", record.ID, start, lastmatch+k))
								}
							}

							c = 0
							start = -1
							flag = true
						} else {
							c++
							if c == 1 { // re-count
								if flag {
									start = i
									nonUniqsNum = 0
									nonUniqs = 0
									lastNonUniqsNum = 0
								}
							}
						}
					} else {
						c++
						if c == 1 { // re-count
							if flag {
								start = i
								nonUniqsNum = 0
								nonUniqs = 0
								lastNonUniqsNum = 0
							}
						}
					}

					if c >= 1 { // at least 1 continuous sites.
						lastmatch = i
						lastNonUniqsNum = nonUniqsNum
					}
				} else { // k-mer not found
					nonUniqs++
					if nonUniqs == 1 {
						nonUniqsNum++
					}
					if nonUniqs <= maxContNonUniqKmers && nonUniqsNum <= maxContNonUniqKmersNum {
						c = 0
						if start > 0 {
							flag = false
						}
					} else {
						if lastNonUniqsNum <= maxContNonUniqKmersNum &&
							start >= 0 && lastmatch-start+k >= minLen {
							if outputFASTA {
								outfh.WriteString(fmt.Sprintf(">%s:%d-%d\n%s\n", record.ID, start+1, lastmatch+k,
									record.Seq.SubSeq(start+1, lastmatch+k).FormatSeq(60)))
							} else {
								outfh.WriteString(fmt.Sprintf("%s\t%d\t%d\n", record.ID, start, lastmatch+k))
							}
						}
						// re-count
						c = 0
						start = -1
						flag = true
					}
				}
				// debug.WriteString(fmt.Sprintln(i, c, start, lastmatch, nonUniqs, nonUniqsNum, lastNonUniqsNum))
			}
			if lastNonUniqsNum <= maxContNonUniqKmersNum+1 &&
				start >= 0 && lastmatch-start+k >= minLen {
				if outputFASTA {
					outfh.WriteString(fmt.Sprintf(">%s:%d-%d\n%s\n", record.ID, start+1, lastmatch+k,
						record.Seq.SubSeq(start+1, lastmatch+k).FormatSeq(60)))
				} else {
					outfh.WriteString(fmt.Sprintf("%s\t%d\t%d\n", record.ID, start, lastmatch+k))
				}
			}
		}
	},
}

func init() {
	RootCmd.AddCommand(uniqsCmd)

	uniqsCmd.Flags().StringP("out-prefix", "o", "-", `out file prefix ("-" for stdout)`)
	uniqsCmd.Flags().BoolP("circular", "", false, "circular genome")
	uniqsCmd.Flags().StringP("genome", "g", "", "genome in (gzipped) fasta file")
	uniqsCmd.Flags().IntP("min-len", "m", 200, "minimum length of subsequence")
	uniqsCmd.Flags().BoolP("allow-muliple-mapped-kmer", "M", false, "allow multiple mapped k-mers")
	uniqsCmd.Flags().BoolP("output-fasta", "a", false, "output fasta format instead of BED3")
	uniqsCmd.Flags().IntP("max-cont-non-uniq-kmers", "x", 0, "max continuous non-unique k-mers")
	uniqsCmd.Flags().IntP("max-num-cont-non-uniq-kmers", "X", 0, "max number of continuous non-unique k-mers")
}
