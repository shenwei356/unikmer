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
	"regexp"
	"runtime"
	"sort"
	"strconv"

	"github.com/pkg/errors"
	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/bio/seqio/fastx"
	"github.com/shenwei356/unikmer"
	"github.com/spf13/cobra"
	"github.com/will-rowe/nthash"
)

var countCmd = &cobra.Command{
	Use:   "count",
	Short: "Count k-mers from FASTA/Q sequences",
	Long: `Count k-mers from FASTA/Q sequences

`,
	Run: func(cmd *cobra.Command, args []string) {
		opt := getOptions(cmd)
		runtime.GOMAXPROCS(opt.NumCPUs)
		seq.ValidateSeq = false

		var err error

		outFile := getFlagString(cmd, "out-prefix")
		circular := getFlagBool(cmd, "circular")
		k := getFlagPositiveInt(cmd, "kmer-len")
		hashed := getFlagBool(cmd, "hash")
		if k > 32 && !hashed {
			hashed = true
			log.Warning("flag -H/--hash is switched on for k > 32")
		}

		canonical := getFlagBool(cmd, "canonical")
		sortKmers := getFlagBool(cmd, "sort")

		if opt.Compact {
			if sortKmers {
				log.Infof("flag -s/--sort overides -c/--compact")
			}
			if hashed {
				log.Infof("flag -H/--hash overides -c/--compact")
			}
		}

		taxid := getFlagUint32(cmd, "taxid")
		var setGlobalTaxid bool

		parseTaxid := getFlagBool(cmd, "parse-taxid")
		parseTaxidRegexp := getFlagString(cmd, "parse-taxid-regexp")

		repeated := getFlagBool(cmd, "repeated")
		moreVerbose := getFlagBool(cmd, "more-verbose")

		if moreVerbose {
			opt.Verbose = true
		}

		var reParseTaxid *regexp.Regexp
		if parseTaxid {
			if taxid > 0 {
				checkError(fmt.Errorf("flag -t/--taxid and -T/--parse-taxid can not given simultaneously"))
			}
			if parseTaxidRegexp == "" {
				checkError(fmt.Errorf("flag -r/--parse-taxid-regexp needed when given flag -T/--parse-taxid"))
			}
			if !regexp.MustCompile(`\(.+\)`).MatchString(parseTaxidRegexp) {
				checkError(fmt.Errorf(`value of -r/--parse-taxid-regexp must contains "(" and ")" to capture taxid`))
			}

			reParseTaxid, err = regexp.Compile(parseTaxidRegexp)
			if err != nil {
				checkError(fmt.Errorf("invalid regular express: %s", parseTaxidRegexp))
			}
		} else if taxid > 0 {
			setGlobalTaxid = true
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

		if !isStdout(outFile) {
			outFile += extDataFile
		}
		outfh, gw, w, err := outStream(outFile, opt.Compress, opt.CompressionLevel)
		checkError(err)
		defer func() {
			outfh.Flush()
			if gw != nil {
				gw.Close()
			}
			w.Close()
		}()

		var mode uint32
		var writer *unikmer.Writer

		if setGlobalTaxid && opt.Verbose {
			log.Infof("set global taxid: %d", taxid)
		}

		if !parseTaxid && !sortKmers {
			if sortKmers {
				mode |= unikmer.UNIK_SORTED
			} else if opt.Compact && !hashed {
				mode |= unikmer.UNIK_COMPACT
			}
			if canonical {
				mode |= unikmer.UNIK_CANONICAL
			}
			if parseTaxid {
				mode |= unikmer.UNIK_INCLUDETAXID
			}
			if hashed {
				mode |= unikmer.UNIK_HASHED
			}
			writer, err = unikmer.NewWriter(outfh, k, mode)
			checkError(err)
			writer.SetMaxTaxid(opt.MaxTaxid)
			if setGlobalTaxid {
				checkError(writer.SetGlobalTaxid(taxid))
			}
		}

		var m map[uint64]struct{}
		var taxondb *unikmer.Taxonomy
		var mt map[uint64]uint32

		// could use bloom filter
		// a key exists means it appear once, value of true means it's appeared more than once.
		var marks map[uint64]bool

		if parseTaxid {
			mt = make(map[uint64]uint32, mapInitSize)
			taxondb = loadTaxonomy(opt, false)
		} else {
			m = make(map[uint64]struct{}, mapInitSize)
		}
		if repeated {
			marks = make(map[uint64]bool, mapInitSize)
		}

		var sequence, kmer, preKmer []byte
		var originalLen, l, end, e int
		var record *fastx.Record
		var fastxReader *fastx.Reader
		var kcode, preKcode unikmer.KmerCode
		var first bool
		var i, j, iters int
		var ok bool
		var n int64
		var founds [][][]byte
		var val uint64
		var lca uint32
		var mark bool
		var nseq int64
		var hasher *nthash.NTHi
		var hash uint64
		for _, file := range files {
			if opt.Verbose {
				log.Infof("reading sequence file: %s", file)
			}
			fastxReader, err = fastx.NewDefaultReader(file)
			checkError(errors.Wrap(err, file))
			for {
				record, err = fastxReader.Read()
				if err != nil {
					if err == io.EOF {
						break
					}
					checkError(errors.Wrap(err, file))
					break
				}

				if parseTaxid {
					founds = reParseTaxid.FindAllSubmatch(record.Name, 1)
					if len(founds) == 0 {
						checkError(fmt.Errorf("failed to parse taxid in header: %s", record.Name))
					}
					val, err = strconv.ParseUint(string(founds[0][1]), 10, 32)
					if err != nil {
						checkError(fmt.Errorf("failed to parse taxid '%s' in header: %s", founds[0][1], record.Name))
					}
					taxid = uint32(val)
				}

				nseq++
				if opt.Verbose && moreVerbose {
					if parseTaxid {
						log.Infof("processing sequence #%d: %s, taxid: %d", nseq, record.ID, taxid)
					} else {
						log.Infof("processing sequence #%d: %s", nseq, record.ID)
					}
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

					// using ntHash
					if hashed {

						hasher, err = nthash.NewHasher(&sequence, uint(k))
						checkError(errors.Wrap(err, file))

						for hash = range hasher.Hash(canonical) {
							if parseTaxid {
								if repeated {
									if mark, ok = marks[hash]; !ok {
										marks[hash] = false
									} else if !mark {
										if lca, ok = mt[hash]; !ok {
											mt[hash] = taxid
										} else {
											mt[hash] = taxondb.LCA(lca, taxid) // update with LCA
										}
										marks[hash] = true
									}

									continue
								}

								if lca, ok = mt[hash]; !ok {
									mt[hash] = taxid
								} else {
									mt[hash] = taxondb.LCA(lca, taxid) // update with LCA
								}
								continue
							}

							if repeated {
								if mark, ok = marks[hash]; !ok {
									marks[hash] = false
								} else if !mark {
									if !sortKmers {
										writer.WriteCode(hash)
										n++
									} else {
										m[hash] = struct{}{}
									}
									marks[hash] = true
								}

								continue
							}

							if _, ok = m[hash]; !ok {
								m[hash] = struct{}{}
								if !sortKmers {
									writer.WriteCode(hash)
									n++
								}
							}
						}

						continue
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

						if canonical {
							kcode = kcode.Canonical()
						}

						if parseTaxid {
							if repeated {
								if mark, ok = marks[kcode.Code]; !ok {
									marks[kcode.Code] = false
								} else if !mark {
									if lca, ok = mt[kcode.Code]; !ok {
										mt[kcode.Code] = taxid
									} else {
										mt[kcode.Code] = taxondb.LCA(lca, taxid) // update with LCA
									}
									marks[kcode.Code] = true
								}

								continue
							}

							if lca, ok = mt[kcode.Code]; !ok {
								mt[kcode.Code] = taxid
							} else {
								mt[kcode.Code] = taxondb.LCA(lca, taxid) // update with LCA
							}
							continue
						}

						if repeated {
							if mark, ok = marks[kcode.Code]; !ok {
								marks[kcode.Code] = false
							} else if !mark {
								if !sortKmers {
									writer.WriteCode(kcode.Code)
									n++
								} else {
									m[kcode.Code] = struct{}{}
								}
								marks[kcode.Code] = true
							}

							continue
						}

						if _, ok = m[kcode.Code]; !ok {
							m[kcode.Code] = struct{}{}
							if !sortKmers {
								writer.WriteCode(kcode.Code)
								n++
							}
						}
					}
				}
			}
		}

		if sortKmers || parseTaxid {
			var mode uint32
			if sortKmers {
				mode |= unikmer.UNIK_SORTED
			} else if opt.Compact && !hashed {
				mode |= unikmer.UNIK_COMPACT
			}
			if canonical {
				mode |= unikmer.UNIK_CANONICAL
			}
			if parseTaxid {
				mode |= unikmer.UNIK_INCLUDETAXID
			}
			if hashed {
				mode |= unikmer.UNIK_HASHED
			}
			writer, err = unikmer.NewWriter(outfh, k, mode)
			checkError(errors.Wrap(err, outFile))
			writer.SetMaxTaxid(opt.MaxTaxid)
			if setGlobalTaxid {
				checkError(writer.SetGlobalTaxid(taxid))
			}

			if parseTaxid {
				n = int64(len(mt))
			} else {
				n = int64(len(m))
			}
			writer.Number = int64(n)

		}

		var code uint64
		if !sortKmers {
			if parseTaxid {
				for code, taxid = range mt {
					writer.WriteCodeWithTaxid(code, taxid)
				}
				n = int64(len(mt))
			}
		} else {
			if parseTaxid {
				codes := make([]uint64, len(mt))

				i := 0
				for code = range mt {
					codes[i] = code
					i++
				}

				if opt.Verbose {
					log.Infof("sorting %d k-mers", len(codes))
				}
				sort.Sort(unikmer.CodeSlice(codes))
				if opt.Verbose {
					log.Infof("done sorting")
				}
				for _, code = range codes {
					writer.WriteCodeWithTaxid(code, mt[code])
				}
				n = int64(len(mt))
			} else {
				codes := make([]uint64, len(m))

				i := 0
				for code = range m {
					codes[i] = code
					i++
				}

				if opt.Verbose {
					log.Infof("sorting %d k-mers", len(codes))
				}
				sort.Sort(unikmer.CodeSlice(codes))
				if opt.Verbose {
					log.Infof("done sorting")
				}

				for _, code = range codes {
					writer.WriteCode(code)
				}
				n = int64(len(m))
			}
		}

		checkError(writer.Flush())
		if opt.Verbose {
			log.Infof("%d unique k-mers saved to %s", n, outFile)
		}
	},
}

func init() {
	RootCmd.AddCommand(countCmd)

	countCmd.Flags().StringP("out-prefix", "o", "-", `out file prefix ("-" for stdout)`)
	countCmd.Flags().IntP("kmer-len", "k", 0, "k-mer length")
	countCmd.Flags().BoolP("circular", "", false, "circular genome")
	countCmd.Flags().BoolP("canonical", "K", false, "only keep the canonical k-mers")
	countCmd.Flags().BoolP("sort", "s", false, helpSort)
	countCmd.Flags().Uint32P("taxid", "t", 0, "global taxid")
	countCmd.Flags().BoolP("parse-taxid", "T", false, `parse taxid from FASTA/Q header`)
	countCmd.Flags().StringP("parse-taxid-regexp", "r", "", `regular expression for passing taxid`)
	countCmd.Flags().BoolP("repeated", "d", false, `only count duplicated k-mers, for removing singleton in FASTQ`)
	countCmd.Flags().BoolP("more-verbose", "V", false, `print extra verbose information`)
	countCmd.Flags().BoolP("hash", "H", false, `save hash of k-mer, automatically on for k>32. This flag overides global flag -c/--compact`)
}
