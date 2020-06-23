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
	"io/ioutil"
	"os"
	"path/filepath"
	"runtime"
	"sort"
	"strings"

	"github.com/shenwei356/unikmer"
	"github.com/shenwei356/util/pathutil"
	"github.com/shenwei356/util/stringutil"
	"github.com/spf13/cobra"
)

// rfilterCmd represents
var rfilterCmd = &cobra.Command{
	Use:   "rfilter",
	Short: "Filter k-mers by taxonomic rank",
	Long: `Filter k-mers by taxonomic rank

Attentions:
  1. flag -L/--lower-than and -H/--higher-than are exclusive, and can be
     used along with -E/--equal-to which values can be different.
  2. a list of pre-ordered ranks is in ~/.unikmer/ranks.txt, you can give
     your list by -r/--rank-file, with one rank per line.
  3. taxids with no rank will be discarded.

`,
	Run: func(cmd *cobra.Command, args []string) {
		opt := getOptions(cmd)
		runtime.GOMAXPROCS(opt.NumCPUs)

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

		rankFile := getFlagString(cmd, "rank-file")
		discardNorank := getFlagBool(cmd, "discard-norank")

		higher := getFlagString(cmd, "higher-than")
		lower := getFlagString(cmd, "lower-than")
		equal := getFlagString(cmd, "equal-to")
		noRank := getFlagString(cmd, "no-rank")

		listOrder := getFlagBool(cmd, "list-order")
		listRanks := getFlagBool(cmd, "list-ranks")

		if higher != "" && lower != "" {
			checkError(fmt.Errorf("-H/--higher-than and -L/--lower-than can't be simultaneous given"))
		}

		rankOrder, err := readRankOrder(opt, rankFile)
		checkError(err)

		if listOrder {
			orders := make([]stringutil.StringCount, 0, len(rankOrder))
			for r, o := range rankOrder {
				orders = append(orders, stringutil.StringCount{Key: r, Count: o})
			}
			sort.Sort(stringutil.ReversedStringCountList{orders})
			for _, order := range orders {
				// fmt.Printf("%d\t%s\n", order.Count, order.Key)
				fmt.Printf("%s\n", order.Key)
			}
			return
		}

		taxondb := loadTaxonomy(opt, true)

		if listRanks {
			orders := make([]stringutil.StringCount, 0, len(taxondb.Ranks))
			var ok bool
			for rank := range taxondb.Ranks {
				if _, ok = rankOrder[rank]; !ok {
					checkError(fmt.Errorf("rank order not defined: %s", rank))
				}
				orders = append(orders, stringutil.StringCount{Key: rank, Count: rankOrder[rank]})
			}
			sort.Sort(stringutil.ReversedStringCountList{orders})
			for _, order := range orders {
				// fmt.Printf("%d\t%s\n", order.Count, order.Key)
				fmt.Printf("%s\n", order.Key)
			}
			return
		}

		filter, err := newRankFilter(taxondb, rankOrder, lower, higher, equal, noRank, discardNorank)
		checkError(err)

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

		var writer *unikmer.Writer

		var infh *bufio.Reader
		var r *os.File
		var reader *unikmer.Reader
		var code uint64
		var taxid uint32
		var k int = -1
		var canonical bool
		var hasTaxid bool
		var flag int
		var nfiles = len(files)
		var n int64
		var rank string
		var pass bool
		for i, file := range files {
			if opt.Verbose {
				log.Infof("processing file (%d/%d): %s", i+1, nfiles, file)
			}

			flag = func() int {
				infh, r, _, err = inStream(file)
				checkError(err)
				defer r.Close()

				reader, err = unikmer.NewReader(infh)
				checkError(err)

				hasTaxid = !opt.IgnoreTaxid && reader.HasTaxidInfo()
				if k == -1 {
					k = reader.K
					canonical = reader.IsCanonical()

					if !hasTaxid {
						checkError(fmt.Errorf(`taxid information not found: %s`, file))
					}

					mode := reader.Flag
					mode |= unikmer.UNIK_INCLUDETAXID
					writer, err = unikmer.NewWriter(outfh, k, mode)
					checkError(err)
					writer.SetMaxTaxid(maxUint32N(reader.GetTaxidBytesLength())) // follow reader
				} else {
					if k != reader.K {
						checkError(fmt.Errorf("K (%d) of binary file '%s' not equal to previous K (%d)", reader.K, file, k))
					}
					if reader.IsCanonical() != canonical {
						checkError(fmt.Errorf(`'canonical' flags not consistent, please check with "unikmer stats"`))
					}
					if !hasTaxid {
						checkError(fmt.Errorf(`taxid information not found: %s`, file))
					}
				}

				for {
					code, taxid, err = reader.ReadCodeWithTaxid()
					if err != nil {
						if err == io.EOF {
							break
						}
						checkError(err)
					}

					rank = taxondb.Rank(taxid)
					if rank == "" {
						continue
					}

					pass, err = filter.rankFilter(rank)
					checkError(err)
					if !pass {
						continue
					}

					n++
					writer.WriteCodeWithTaxid(code, taxid) // not need to check err
					// fmt.Printf("%d\t%s\n", taxid, rank)
				}

				return flagContinue
			}()

			if flag == flagReturn {
				return
			} else if flag == flagBreak {
				break
			}
		}

		checkError(writer.Flush())
		if opt.Verbose {
			log.Infof("%d k-mers saved to %s", n, outFile)
		}
	},
}

func init() {
	RootCmd.AddCommand(rfilterCmd)

	rfilterCmd.Flags().StringP("out-prefix", "o", "-", `out file prefix ("-" for stdout)`)
	rfilterCmd.Flags().StringP("rank-file", "r", "", "user-defined ordered taxonomic ranks")
	rfilterCmd.Flags().BoolP("list-order", "", false, "list defined ranks in order")
	rfilterCmd.Flags().BoolP("list-ranks", "", false, "list ordered ranks in taxonomy database")

	rfilterCmd.Flags().BoolP("discard-norank", "n", false, `discard "no rank", defined by --no-rank`)
	rfilterCmd.Flags().StringP("no-rank", "", "no rank", `value of "no rank"`)

	rfilterCmd.Flags().StringP("lower-than", "L", "", "output ranks lower than a rank, exclusive with --higher-than")
	rfilterCmd.Flags().StringP("higher-than", "H", "", "output ranks higher than a rank, exclusive with --lower-than")
	rfilterCmd.Flags().StringP("equal-to", "E", "", "output ranks equal to a rank")
}

type rankFilter struct {
	db        *unikmer.Taxonomy
	rankOrder map[string]int

	lower  string
	higher string
	equal  string

	oLower  int
	oHigher int
	oEqual  int

	limitLower  bool
	limitHigher bool
	limitEqual  bool

	noRank       string
	discardNoank bool
}

func newRankFilter(db *unikmer.Taxonomy, rankOrder map[string]int, lower, higher, equal, noRank string, discardNorank bool) (*rankFilter, error) {
	if lower != "" && higher != "" {
		return nil, fmt.Errorf("higher and lower can't be simultaneous given")
	}
	f := &rankFilter{
		db:           db,
		rankOrder:    rankOrder,
		lower:        lower,
		higher:       higher,
		equal:        equal,
		noRank:       noRank,
		discardNoank: discardNorank,
	}
	var err error
	if lower != "" {
		f.oLower, err = getRankOrder(db, rankOrder, lower)
		if err != nil {
			return nil, err
		}
		f.limitLower = true
	}
	if higher != "" {
		f.oHigher, err = getRankOrder(db, rankOrder, higher)
		if err != nil {
			return nil, err
		}
		f.limitHigher = true
	}
	if equal != "" {
		f.oEqual, err = getRankOrder(db, rankOrder, equal)
		if err != nil {
			return nil, err
		}
		f.limitEqual = true
	}
	if noRank != "" {
		_, err = getRankOrder(db, rankOrder, noRank)
		if err != nil {
			return nil, err
		}
	}
	return f, nil
}

func getRankOrder(db *unikmer.Taxonomy, rankOrder map[string]int, rank string) (int, error) {
	var ok bool
	if _, ok = rankOrder[rank]; !ok {
		return -1, fmt.Errorf("rank order not defined in rank file: %s", rank)
	}
	if _, ok = db.Ranks[rank]; !ok {
		return -1, fmt.Errorf("rank order not found in taxonomy database: %s", rank)
	}

	return rankOrder[rank], nil
}

func (f *rankFilter) rankFilter(rank string) (bool, error) {
	pass := false
	order, err := getRankOrder(f.db, f.rankOrder, rank)
	if err != nil {
		return false, fmt.Errorf("query rank: %s", err)
	}

	if f.limitEqual {
		if f.oEqual == order {
			pass = true
		} else if f.limitLower {
			pass = order < f.oLower
		} else if f.limitHigher {
			pass = order > f.oHigher
		} else {
			pass = false
		}
	} else if f.limitLower {
		pass = order < f.oLower
	} else if f.limitHigher {
		pass = order > f.oHigher
	} else {
		pass = true // no any filter
	}

	if pass && f.discardNoank && rank == f.noRank {
		pass = false
	}

	return pass, nil
}

func readRankOrderFromFile(file string) (map[string]int, error) {
	fh, err := os.Open(file)
	if err != nil {
		return nil, fmt.Errorf("read rank order list from '%s': %s", file, err)
	}

	ranks := make([]string, 0, 100)

	scanner := bufio.NewScanner(fh)
	var rank string
	for scanner.Scan() {
		rank = strings.TrimSpace(scanner.Text())
		if rank == "" || rank[0] == '#' {
			continue
		}

		ranks = append(ranks, rank)
	}
	if err = scanner.Err(); err != nil {
		return nil, fmt.Errorf("read rank order list from '%s': %s", file, err)
	}

	if len(ranks) == 0 {
		return nil, fmt.Errorf("no ranks found in file: %s", file)
	}
	rankOrder := make(map[string]int, len(ranks))
	order := 0
	var ok bool
	for i := len(ranks) - 1; i >= 0; i-- {
		if _, ok = rankOrder[ranks[i]]; ok {
			return nil, fmt.Errorf("duplicated rank: %s", ranks[i])
		}
		rankOrder[ranks[i]] = order
		order++
	}
	return rankOrder, nil
}

func readRankOrder(opt *Options, rankFile string) (map[string]int, error) {
	if rankFile != "" {
		if opt.Verbose {
			log.Infof("read rank order from: %s", rankFile)
		}
		return readRankOrderFromFile(rankFile)
	}

	defaultRankFile := filepath.Join(opt.DataDir, defaultRanksFile)
	existed, err := pathutil.Exists(defaultRankFile)
	if err != nil {
		return nil, fmt.Errorf("check default rank file: %s", defaultRankFile)
	}
	if !existed {
		if opt.Verbose {
			log.Infof("write default rank order to: %s", defaultRankFile)
		}
		err = writeDefaltRankOrderFile(defaultRankFile)
		if err != nil {
			return nil, fmt.Errorf("write default rank file: %s", defaultRankFile)
		}
	}

	if opt.Verbose {
		log.Infof("read rank order from: %s", defaultRankFile)
	}
	return readRankOrderFromFile(defaultRankFile)
}

func writeDefaltRankOrderFile(file string) error {
	return ioutil.WriteFile(file, []byte(defaultRanksText), 0644)
}

const defaultRanksFile = "ranks.txt"
const defaultRanksText = `
# Ref: https://en.wikipedia.org/wiki/Taxonomic_rank
hyperkingdom
superkingdom
kingdom
subkingdom
infrakingdom
parvkingdom

superphylum
phylum
subphylum
infraphylum
microphylum

superclass
class
subclass
infraclass
parvclass

superdivision
division
subdivision
infradivision

superlegion
legion
sublegion
infralegion

supercohort
cohort
subcohort
infracohort

gigaorder
magnorder
grandorder
mirorder
superorder
# series
order
# parvorder
nanorder
hypoorder
minorder
suborder
infraorder
parvorder

# section
# subsection

gigafamily
megafamily
grandfamily
hyperfamily
superfamily
epifamily
# series
group
family
subfamily
infrafamily

supertribe
tribe
subtribe
infratribe

genus
subgenus
section
subsection
series
subseries

species group
species subgroup

superspecies
species
subspecies
varietas
variety
subvarietas
subvariety
forma
form
subforma
subform

no rank

`
