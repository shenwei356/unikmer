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
	"bufio"
	"fmt"
	"io"
	"os"
	"regexp"
	"strings"

	"github.com/pkg/errors"
	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/bio/seqio/fastx"
	"github.com/shenwei356/bio/sketches"
	"github.com/shenwei356/unik/v5"

	"github.com/spf13/cobra"
)

var mapCmd = &cobra.Command{
	Use:     "map",
	Aliases: []string{"uniqs"},
	Short:   "Mapping k-mers back to the genome and extract successive regions/subsequences",
	Long: `Mapping k-mers back to the genome and extract successive regions/subsequences

Attention:
  0. By default, only unique-mapped k-mers are considered.
     You can use -M/--allow-multiple-mapped-kmerss to allow mutiple-mapped k-mers.
  1. The 'canonical/scaled/hashed' flags of all files should be consistent.
  2. Default output is in BED3 format, with left-closed and right-open
     0-based interval.
  3. When using flag --circular, end position of subsequences that 
     crossing genome sequence end would be greater than sequence length.

`,
	Run: func(cmd *cobra.Command, args []string) {
		opt := getOptions(cmd)
		seq.ValidateSeq = false

		var err error

		reSeqNameStrs := getFlagStringSlice(cmd, "seq-name-filter")
		reSeqNames := make([]*regexp.Regexp, 0, len(reSeqNameStrs))
		for _, kw := range reSeqNameStrs {
			if !reIgnoreCase.MatchString(kw) {
				kw = reIgnoreCaseStr + kw
			}
			re, err := regexp.Compile(kw)
			if err != nil {
				checkError(errors.Wrapf(err, "failed to parse regular expression for matching sequence header: %s", kw))
			}
			reSeqNames = append(reSeqNames, re)
		}
		filterNames := len(reSeqNames) > 0

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

		checkFileSuffix(opt, extDataFile, files...)

		outFile := getFlagString(cmd, "out-prefix")

		genomes := getFlagStringSlice(cmd, "genome")
		if len(genomes) == 0 {
			checkError(fmt.Errorf("flag -g/--genome needed"))
		}

		minLen := getFlagPositiveInt(cmd, "min-len")
		mMapped := getFlagBool(cmd, "allow-multiple-mapped-kmers")
		outputFASTA := getFlagBool(cmd, "output-fasta")
		maxGapSize := getFlagNonNegativeInt(cmd, "max-gap-size")
		maxGapNum := getFlagNonNegativeInt(cmd, "max-gap-num")
		seqsAsOneGenome := getFlagBool(cmd, "seqs-in-a-file-as-one-genome")
		circular := getFlagBool(cmd, "circular")

		if seqsAsOneGenome && mMapped {
			checkError(fmt.Errorf("flag -M/--allow-multiple-mapped-kmers and -W/--seqs-in-a-file-as-one-genome are not compatible"))
		}

		if maxGapNum > 0 && maxGapSize == 0 {
			log.Warningf("-X/--max-gap-num %d is ignored becaue value of -x/--max-gap-size is 0", maxGapNum)
		}
		if maxGapSize > 0 && maxGapNum == 0 {
			checkError(fmt.Errorf("value of -X/--max-gap-num should be > 0 when value of -x/--max-gap-size is > 0"))
		}

		// all kmers in .unik files
		m := make(map[uint64]struct{}, mapInitSize)

		// -----------------------------------------------------------------------

		var k int = -1
		var canonical bool
		var infh *bufio.Reader
		var r *os.File
		var reader0 *unik.Reader
		var hashed bool
		var code uint64
		var nfiles = len(files)
		for i, file := range files {
			if opt.Verbose {
				log.Infof("reading file (%d/%d): %s", i+1, nfiles, file)
			}
			func() {
				infh, r, _, err = inStream(file)
				checkError(err)
				defer r.Close()

				reader, err := unik.NewReader(infh)
				checkError(errors.Wrap(err, file))

				if k == -1 {
					reader0 = reader
					k = reader.K
					hashed = reader.IsHashed()
					canonical = reader.IsCanonical()
					if !canonical {
						checkError(fmt.Errorf("%s: 'canonical' flag is needed", file))
					}
				} else {
					checkCompatibility(reader0, reader, file)
				}

				for {
					code, _, err = reader.ReadCodeWithTaxid()
					if err != nil {
						if err == io.EOF {
							break
						}
						checkError(errors.Wrap(err, file))
					}

					m[code] = struct{}{}
				}
			}()
		}

		if opt.Verbose {
			log.Infof("%d k-mers loaded", len(m))
		}

		// -----------------------------------------------------------------------
		var m2 map[int]map[uint64]bool // genome-id -> kmer -> mutiple-mapped
		var _m2 map[uint64]bool

		var record *fastx.Record
		var fastxReader *fastx.Reader
		var iter *sketches.Iterator
		var i int
		var ok bool
		var multipleMapped bool
		var ignoreSeq bool
		var re *regexp.Regexp

		if !mMapped {
			m2 = make(map[int]map[uint64]bool, 8)
			var nKmers uint64
			var genomeIdx int
			for _, genomeFile := range genomes {
				if opt.Verbose {
					log.Infof("pre-reading genome file: %s", genomeFile)
				}

				fastxReader, err = fastx.NewDefaultReader(genomeFile)
				checkError(errors.Wrap(err, genomeFile))
				for {
					record, err = fastxReader.Read()
					if err != nil {
						if err == io.EOF {
							break
						}
						checkError(errors.Wrap(err, genomeFile))
						break
					}

					if filterNames {
						ignoreSeq = false
						for _, re = range reSeqNames {
							if re.Match(record.Name) {
								ignoreSeq = true
								break
							}
						}
						if ignoreSeq {
							continue
						}
					}

					if len(record.Seq.Seq) < k {
						continue
					}

					if hashed {
						iter, err = sketches.NewHashIterator(record.Seq, k, true, circular)
					} else {
						iter, err = sketches.NewKmerIterator(record.Seq, k, true, circular)
					}
					if err != nil {
						if err == sketches.ErrShortSeq {
							if opt.Verbose {
								log.Infof("ignore short seq in file '%s': %s", genomeFile, record.Name)
							}
							continue
						} else {
							checkError(errors.Wrapf(err, "file: %s, seq: %s", genomeFile, record.Name))
						}
					}

					if _m2, ok = m2[genomeIdx]; !ok {
						_m2 = make(map[uint64]bool, mapInitSize)
						m2[genomeIdx] = _m2
					}

					for {
						code, ok, err = iter.Next()
						if hashed && err != nil {
							checkError(errors.Wrapf(err, "%s: %s", record.Name, record.Seq.Seq[iter.Index():iter.Index()+k]))
						}
						if !ok {
							break
						}

						if multipleMapped, ok = _m2[code]; !ok {
							nKmers++
							_m2[code] = false
						} else if !multipleMapped {
							_m2[code] = true
						}
					}

					if !seqsAsOneGenome {
						genomeIdx++
					}
				}
			}

			if opt.Verbose {
				log.Infof("%d k-mers loaded from %d genomes", nKmers, len(m2))
			}

			for genomeIdx, _m2 = range m2 {
				for code, multipleMapped = range _m2 {
					if !multipleMapped {
						delete(m2[genomeIdx], code)
					}
				}
			}
			n := 0
			for _, _m2 = range m2 {
				n += len(_m2)
			}
			if opt.Verbose {
				log.Infof("%d k-mers in genomes are multiple mapped", n)
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

		var genomeIdx int
		for _, genomeFile := range genomes {
			var c, start, nonUniqs, nonUniqsNum, lastNonUniqsNum, lastmatch int
			var length0 int // origninal length of sequence
			var flag bool = true
			if opt.Verbose {
				log.Infof("reading genome file: %s", genomeFile)
			}
			fastxReader, err = fastx.NewDefaultReader(genomeFile)
			checkError(errors.Wrap(err, genomeFile))
			for {
				record, err = fastxReader.Read()
				if err != nil {
					if err == io.EOF {
						break
					}
					checkError(errors.Wrap(err, genomeFile))
					break
				}

				if filterNames {
					ignoreSeq = false
					for _, re = range reSeqNames {
						if re.Match(record.Name) {
							ignoreSeq = true
							break
						}
					}
					if ignoreSeq {
						continue
					}
				}

				if len(record.Seq.Seq) < k {
					continue
				}

				length0 = len(record.Seq.Seq)

				if circular { // concat two copies of sequence
					record.Seq.Seq = append(record.Seq.Seq, record.Seq.Seq...)
				}

				if opt.Verbose {
					log.Infof("processinig sequence: %s", record.ID)
				}

				c = 0
				start = -1
				nonUniqs = 0
				nonUniqsNum = 0

				if hashed {
					iter, err = sketches.NewHashIterator(record.Seq, k, true, false)
				} else {
					iter, err = sketches.NewKmerIterator(record.Seq, k, true, false)
				}
				checkError(errors.Wrapf(err, "seq: %s", record.Name))

				if !mMapped {
					_m2 = m2[genomeIdx]
				}

				for {
					code, ok, err = iter.Next()
					if !hashed && err != nil {
						checkError(errors.Wrapf(err, "%s: %s", record.Name, record.Seq.Seq[iter.Index():iter.Index()+k]))
					}
					if !ok {
						break
					}

					i = iter.Index()

					if _, ok = m[code]; ok {
						nonUniqs = 0
						if !mMapped {
							if multipleMapped, ok = _m2[code]; ok && multipleMapped {
								if lastNonUniqsNum <= maxGapNum &&
									start >= 0 && lastmatch-start+k >= minLen {

									// subsequence longer than original sequence
									if circular && lastmatch-start+k > length0 {
										lastmatch = length0 - k + start
									}

									if outputFASTA {
										fmt.Fprintf(outfh, ">%s:%d-%d\n%s\n", record.ID, start+1, lastmatch+k,
											record.Seq.SubSeq(start+1, lastmatch+k).FormatSeq(60))
									} else {
										fmt.Fprintf(outfh, "%s\t%d\t%d\n", record.ID, start, lastmatch+k)
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

										// 2nd clone of seq
										if circular && start >= length0 {
											break
										}
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

									// 2nd clone of seq
									if circular && start >= length0 {
										break
									}
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
						if nonUniqs <= maxGapSize && nonUniqsNum <= maxGapNum {
							c = 0
							if start > 0 {
								flag = false
							}
						} else {
							if lastNonUniqsNum <= maxGapNum &&
								start >= 0 && lastmatch-start+k >= minLen {

								// subsequence longer than original sequence
								if circular && lastmatch-start+k > length0 {
									lastmatch = length0 - k + start
								}

								if outputFASTA {
									fmt.Fprintf(outfh, ">%s:%d-%d\n%s\n", record.ID, start+1, lastmatch+k,
										record.Seq.SubSeq(start+1, lastmatch+k).FormatSeq(60))
								} else {
									fmt.Fprintf(outfh, "%s\t%d\t%d\n", record.ID, start, lastmatch+k)
								}
							}
							// re-count
							c = 0
							start = -1
							flag = true
						}
					}

					if !seqsAsOneGenome {
						genomeIdx++
					}

					// debug.WriteString(fmt.Sprintln(i, c, start, lastmatch, nonUniqs, nonUniqsNum, lastNonUniqsNum))
				}
				if lastNonUniqsNum <= maxGapNum+1 &&
					start >= 0 && lastmatch-start+k >= minLen {

					// subsequence longer than original sequence
					if circular && lastmatch-start+k > length0 {
						lastmatch = length0 - k + start
					}

					if outputFASTA {
						fmt.Fprintf(outfh, ">%s:%d-%d\n%s\n", record.ID, start+1, lastmatch+k,
							record.Seq.SubSeq(start+1, lastmatch+k).FormatSeq(60))
					} else {
						fmt.Fprintf(outfh, "%s\t%d\t%d\n", record.ID, start, lastmatch+k)
					}
				}
			}
		}
	},
}

func init() {
	RootCmd.AddCommand(mapCmd)

	mapCmd.Flags().StringSliceP("seq-name-filter", "B", []string{}, `list of regular expressions for filtering out sequences by header/name, case ignored`)

	mapCmd.Flags().StringP("out-prefix", "o", "-", `out file prefix ("-" for stdout)`)
	mapCmd.Flags().StringSliceP("genome", "g", []string{}, "genomes in (gzipped) fasta file(s)")
	mapCmd.Flags().IntP("min-len", "m", 200, "minimum length of subsequence")
	mapCmd.Flags().BoolP("allow-multiple-mapped-kmers", "M", false, "allow multiple mapped k-mers")
	mapCmd.Flags().BoolP("seqs-in-a-file-as-one-genome", "W", false, "treat seqs in a genome file as one genome")
	mapCmd.Flags().BoolP("output-fasta", "a", false, "output fasta format instead of BED3")

	mapCmd.Flags().IntP("max-gap-size", "x", 0, "max gap size (the number of consecutive unmapped k-mers)")
	mapCmd.Flags().IntP("max-gap-num", "X", 0, "max number of gaps (consecutive unmapped k-mers)")
	mapCmd.Flags().BoolP("circular", "", false, `circular genome. type "unikmer uniqs -h" for details`)
}
