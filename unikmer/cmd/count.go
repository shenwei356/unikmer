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
	"io"
	"regexp"
	"strconv"
	"strings"

	"github.com/pkg/errors"
	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/bio/seqio/fastx"
	"github.com/shenwei356/bio/sketches"
	"github.com/shenwei356/bio/taxdump"
	"github.com/shenwei356/unik/v5"

	"github.com/spf13/cobra"
	"github.com/twotwotwo/sorts/sortutil"
)

var countCmd = &cobra.Command{
	Use:   "count",
	Short: "Generate k-mers (sketch) from FASTA/Q sequences",
	Long: `Generate k-mers (sketch) from FASTA/Q sequences

K-mer:
  1. K-mer code (k<=32)
  2. Hashed k-mer (ntHash, k<=64)

K-mer sketches:
  1. Scaled MinHash
  2. Minimizer
  3. Closed Syncmer

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

		outFile := getFlagString(cmd, "out-prefix")
		k := getFlagPositiveInt(cmd, "kmer-len")
		canonical := getFlagBool(cmd, "canonical")

		hashed := getFlagBool(cmd, "hash")
		if k > 32 && !hashed {
			hashed = true
			log.Warning("flag -H/--hash is switched on for k > 32")
		}
		if hashed && k > 64 {
			checkError(fmt.Errorf("k-mer size (%d) should be <=64", k))
		}

		scale := getFlagPositiveInt(cmd, "scale")
		if scale > 1<<31-1 {
			checkError(fmt.Errorf("value of flag --scale is too big"))
		}
		scaled := scale > 1
		if scaled && !hashed {
			hashed = true
			log.Warning("flag -H/--hash is switched on for scale > 1")
		}
		maxHash := uint64(float64(^uint64(0)) / float64(scale))

		minimizerW := getFlagNonNegativeInt(cmd, "minimizer-w")
		if minimizerW > 1<<31-1 {
			checkError(fmt.Errorf("value of flag --minimizer-w is too big"))
		}
		minimizer := minimizerW > 0
		if minimizer {
			if !hashed {
				hashed = true
				log.Warning("flag -H/--hash is switched on for minimizer-w > 1")
			}
			if !canonical {
				hashed = true
				log.Warning("flag -K/--canonical is switched on for minimizer-w > 1")
			}
		}

		syncmerS := getFlagNonNegativeInt(cmd, "syncmer-s")
		if syncmerS > 1<<31-1 {
			checkError(fmt.Errorf("value of flag --syncmer-s is too big"))
		}
		syncmer := syncmerS > 0
		if syncmer {
			if !hashed {
				hashed = true
				log.Warning("flag -H/--hash is switched on for syncmer-s > 1")
			}
			if !canonical {
				hashed = true
				log.Warning("flag -K/--canonical is switched on for syncmer-s > 1")
			}
		}
		if minimizer && syncmer {
			checkError(fmt.Errorf("flag --minimizer-w and --syncmer-s can not be given simultaneously"))
		}

		sortKmers := getFlagBool(cmd, "sort")
		circular := getFlagBool(cmd, "circular")

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
		unique := getFlagBool(cmd, "unique")

		if unique && repeated {
			checkError(fmt.Errorf("flag -u/--unique and -d/--repeated are not compatible"))
		}

		linear := getFlagBool(cmd, "linear")
		if linear && (repeated || unique) {
			log.Warningf("flag -d/repeated and -u/--unique are ignored when -l/--linear given")
		}
		if linear && sortKmers {
			checkError(fmt.Errorf("flag -l/--linear and -s/--sort are not compatible"))
		}

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
		files := getFileListFromArgsAndFile(cmd, args, true, "infile-list", !opt.SkipFileCheck)
		if opt.Verbose {
			if len(files) == 1 && isStdin(files[0]) {
				log.Info("no files given, reading from stdin")
			} else {
				log.Infof("%d input file(s) given", len(files))
			}
		}

		if !isStdout(outFile) && !strings.HasSuffix(outFile, extDataFile) {
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

		if setGlobalTaxid && opt.Verbose {
			log.Infof("set global taxid: %d", taxid)
		}

		var writer *unik.Writer
		var mode uint32
		var n uint64

		var m map[uint64]struct{}
		var taxondb *taxdump.Taxonomy
		var mt map[uint64]uint32

		// could use bloom filter
		// a key exists means it appear once, value of true means it's appeared more than once.
		var marks map[uint64]bool

		if linear {
			if opt.Compact && !hashed {
				mode |= unik.UnikCompact
			}
			if canonical {
				mode |= unik.UnikCanonical
			}
			if parseTaxid {
				mode |= unik.UnikIncludeTaxID
			}
			if hashed {
				mode |= unik.UnikHashed
			}
			writer, err = unik.NewWriter(outfh, k, mode)
			checkError(errors.Wrap(err, outFile))
			writer.SetMaxTaxid(opt.MaxTaxid)
			if setGlobalTaxid {
				checkError(writer.SetGlobalTaxid(taxid))
			}
			if scaled {
				writer.SetScale(uint32(scale))
			}

			n = 0
		} else {
			if parseTaxid {
				mt = make(map[uint64]uint32, mapInitSize)
				taxondb = loadTaxonomy(opt, false)
			} else if !(repeated || unique) {
				m = make(map[uint64]struct{}, mapInitSize)
			}
			if repeated || unique {
				marks = make(map[uint64]bool, mapInitSize)
			}
		}

		var record *fastx.Record
		var fastxReader *fastx.Reader
		var ok bool
		var founds [][][]byte
		var val uint64
		var lca uint32
		var mark bool
		var nseq int64
		var code uint64
		var iter *sketches.Iterator
		var sketch *sketches.Sketch
		var ignoreSeq bool
		var re *regexp.Regexp

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

				if syncmer {
					sketch, err = sketches.NewSyncmerSketch(record.Seq, k, syncmerS, circular)
				} else if minimizer {
					sketch, err = sketches.NewMinimizerSketch(record.Seq, k, minimizerW, circular)
				} else if hashed {
					iter, err = sketches.NewHashIterator(record.Seq, k, canonical, circular)
				} else {
					iter, err = sketches.NewKmerIterator(record.Seq, k, canonical, circular)
				}
				if err != nil {
					if err == sketches.ErrShortSeq {
						if opt.Verbose && moreVerbose {
							log.Infof("ignore short seq: %s", record.Name)
						}
						continue
					} else {
						checkError(errors.Wrapf(err, "seq: %s", record.Name))
					}
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

				for {
					if syncmer {
						code, ok = sketch.NextSyncmer()
					} else if minimizer {
						code, ok = sketch.NextMinimizer()
					} else if hashed {
						code, ok = iter.NextHash()
					} else {
						code, ok, err = iter.NextKmer()
						if err != nil {
							checkError(errors.Wrapf(err, "seq: %s", record.Name))
						}
					}

					if !ok {
						break
					}

					if scaled && code > maxHash {
						continue
					}

					if parseTaxid {
						if repeated {
							if mark, ok = marks[code]; !ok {
								mt[code] = taxid
								marks[code] = false
							} else {
								if lca, ok = mt[code]; !ok {
									mt[code] = taxid
								} else {
									mt[code] = taxondb.LCA(lca, taxid) // update with LCA
								}
								if !mark {
									marks[code] = true
								}
							}

							continue
						} else if unique {
							if mark, ok = marks[code]; !ok {
								mt[code] = taxid // though added here, but can't ensure it's uniuqe.
								marks[code] = false
							} else if !mark {
								marks[code] = true
							}

							continue
						}

						if lca, ok = mt[code]; !ok {
							mt[code] = taxid
						} else {
							mt[code] = taxondb.LCA(lca, taxid) // update with LCA
						}
						continue
					}

					if linear {
						if parseTaxid {
							writer.WriteCodeWithTaxid(code, taxid)
						} else {
							writer.WriteCode(code)
						}
						n++

						continue
					}

					if repeated || unique {
						if mark, ok = marks[code]; !ok {
							marks[code] = false
						} else if !mark {
							marks[code] = true
						}

						continue
					}

					if _, ok = m[code]; !ok {
						m[code] = struct{}{}
					}
				}
			}
		}

		if linear {
			checkError(writer.Flush())
			if opt.Verbose {
				log.Infof("%d unique k-mers saved to %s", n, outFile)
			}
			return
		}

		if sortKmers {
			mode |= unik.UnikSorted
		} else if opt.Compact && !hashed {
			mode |= unik.UnikCompact
		}
		if canonical {
			mode |= unik.UnikCanonical
		}
		if parseTaxid {
			mode |= unik.UnikIncludeTaxID
		}
		if hashed {
			mode |= unik.UnikHashed
		}
		writer, err = unik.NewWriter(outfh, k, mode)
		checkError(errors.Wrap(err, outFile))
		writer.SetMaxTaxid(opt.MaxTaxid)
		if setGlobalTaxid {
			checkError(writer.SetGlobalTaxid(taxid))
		}
		if scaled {
			writer.SetScale(uint32(scale))
		}

		n = 0

		if repeated {
			for _, mark = range marks {
				if mark {
					n++
				}
			}
		} else if unique {
			for _, mark = range marks {
				if !mark {
					n++
				}
			}
		} else if parseTaxid {
			n = uint64(len(mt))
		} else {
			n = uint64(len(m))
		}

		writer.Number = uint64(n)

		if !sortKmers {
			if parseTaxid {
				if repeated {
					for code, mark = range marks {
						if mark {
							writer.WriteCodeWithTaxid(code, mt[code])
						}
					}
				} else if unique {
					for code, mark = range marks {
						if !mark {
							writer.WriteCodeWithTaxid(code, mt[code])
						}
					}
				} else {
					for code, taxid = range mt {
						writer.WriteCodeWithTaxid(code, taxid)
					}
				}
			} else if repeated {
				for code, mark = range marks {
					if mark {
						writer.WriteCode(code)
					}
				}
			} else if unique {
				for code, mark = range marks {
					if !mark {
						writer.WriteCode(code)
					}
				}
			} else {
				for code = range m {
					writer.WriteCode(code)
				}
			}
		} else {
			codes := make([]uint64, n)

			i := 0
			if parseTaxid {
				if repeated {
					for code, mark = range marks {
						if mark {
							codes[i] = code
							i++
						}
					}
				} else if unique {
					for code, mark = range marks {
						if !mark {
							codes[i] = code
							i++
						}
					}
				} else {
					for code = range mt {
						codes[i] = code
						i++
					}
				}
			} else if repeated {
				for code, mark = range marks {
					if mark {
						codes[i] = code
						i++
					}
				}
			} else if unique {
				for code, mark = range marks {
					if !mark {
						codes[i] = code
						i++
					}
				}
			} else {
				for code = range m {
					codes[i] = code
					i++
				}
			}

			if opt.Verbose {
				log.Infof("sorting %d k-mers", len(codes))
			}
			// sort.Sort(sketches.CodeSlice(codes))
			sortutil.Uint64s(codes)
			if opt.Verbose {
				log.Infof("done sorting")
			}

			if parseTaxid {
				for _, code = range codes {
					writer.WriteCodeWithTaxid(code, mt[code])
				}
			} else {
				for _, code = range codes {
					writer.WriteCode(code)
				}
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

	countCmd.Flags().StringSliceP("seq-name-filter", "B", []string{}, `list of regular expressions for filtering out sequences by header/name, case ignored.`)

	countCmd.Flags().StringP("out-prefix", "o", "-", `out file prefix ("-" for stdout)`)
	countCmd.Flags().IntP("kmer-len", "k", 0, "k-mer length")
	countCmd.Flags().BoolP("canonical", "K", false, "only keep the canonical k-mers")
	countCmd.Flags().BoolP("sort", "s", false, helpSort)
	countCmd.Flags().Uint32P("taxid", "t", 0, "global taxid")
	countCmd.Flags().BoolP("parse-taxid", "T", false, `parse taxid from FASTA/Q header`)
	countCmd.Flags().StringP("parse-taxid-regexp", "r", "", `regular expression for passing taxid`)
	countCmd.Flags().BoolP("repeated", "d", false, `only count duplicate k-mers, for removing singleton in FASTQ`)
	countCmd.Flags().BoolP("unique", "u", false, `only count unique k-mers, which are not duplicate`)
	countCmd.Flags().BoolP("more-verbose", "V", false, `print extra verbose information`)
	countCmd.Flags().BoolP("hash", "H", false, `save hash of k-mer, automatically on for k>32. This flag overides global flag -c/--compact`)
	countCmd.Flags().BoolP("circular", "", false, "circular genome")

	countCmd.Flags().IntP("scale", "D", 1, `scale/down-sample factor`)
	countCmd.Flags().IntP("minimizer-w", "W", 0, `minimizer window size`)
	countCmd.Flags().IntP("syncmer-s", "S", 0, `closed syncmer length`)

	countCmd.Flags().BoolP("linear", "l", false, `output k-mers in linear order, duplicate k-mers are not removed`)

	countCmd.SetUsageTemplate(usageTemplate("-K -k <k> -u -s [-t <taxid>] <seq files> -o <out prefix>"))

}

var reIgnoreCaseStr = "(?i)"
var reIgnoreCase = regexp.MustCompile(`\(\?i\)`)
