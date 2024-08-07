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
	"path/filepath"
	"strings"

	"github.com/pkg/errors"

	"github.com/shenwei356/bio/taxdump"
	"github.com/shenwei356/unik/v5"
	"github.com/shenwei356/util/pathutil"
	"github.com/shenwei356/util/stringutil"
	"github.com/spf13/cobra"
	"github.com/twotwotwo/sorts"
)

// rfilterCmd represents
var rfilterCmd = &cobra.Command{
	Use:   "rfilter",
	Short: "Filter k-mers by taxonomic rank",
	Long: `Filter k-mers by taxonomic rank

Attentions:
  1. Flag -L/--lower-than and -H/--higher-than are exclusive, and can be
     used along with -E/--equal-to which values can be different.
  2. A list of pre-ordered ranks is in ~/.unikmer/ranks.txt, you can use
     your list by -r/--rank-file, the format specification is below.
  3. All ranks in taxonomy database should be defined in rank file.
  4. Ranks can be removed with black list via -B/--black-list.
  5. TaxIds with no rank can be optionally discarded by -N/--discard-noranks.
  6. But when filtering with -L/--lower-than, you can use
    -n/--save-predictable-norank to save some special ranks without order,
    where rank of the closest higher node is still lower than rank cutoff.

Rank file:
  1. Blank lines or lines starting with "#" are ignored.
  2. Ranks are in decending order and case ignored.
  3. Ranks with same order should be in one line separated with comma (",", no space).
  4. Ranks without order should be assigned a prefix symbol "!" for each rank.

`,
	Run: func(cmd *cobra.Command, args []string) {
		opt := getOptions(cmd)

		var err error

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

		checkFileSuffix(opt, extDataFile, files...)

		outFile := getFlagString(cmd, "out-prefix")

		rankFile := getFlagString(cmd, "rank-file")

		discardNoRank := getFlagBool(cmd, "discard-noranks")
		saveNorank := getFlagBool(cmd, "save-predictable-norank")

		blackListRanks := getFlagStringSlice(cmd, "black-list")

		rootTaxid := getFlagUint32(cmd, "root-taxid")
		discardRoot := getFlagBool(cmd, "discard-root")

		higher := strings.ToLower(getFlagString(cmd, "higher-than"))
		lower := strings.ToLower(getFlagString(cmd, "lower-than"))
		equalsS := getFlagStringSlice(cmd, "equal-to")
		equals := make([]string, 0, len(equalsS))
		for _, val := range equalsS {
			equals = append(equals, strings.ToLower(val))
		}

		listOrder := getFlagBool(cmd, "list-order")
		listRanks := getFlagBool(cmd, "list-ranks")

		if higher != "" && lower != "" {
			checkError(fmt.Errorf("-H/--higher-than and -L/--lower-than can't be simultaneous given"))
		}

		if saveNorank {
			discardNoRank = true
			if !discardNoRank {
				log.Infof("flag -N/--discard-noranks is switched on when using -n/--save-predictable-norank")
			}

			if lower == "" {
				checkError(fmt.Errorf("flag -n/--save-predictable-norank only works along with -L/--lower-than"))
			}
		}

		rankOrder, noRanks, err := readRankOrder(opt, rankFile)
		checkError(errors.Wrap(err, rankFile))

		noRanksList := make([]string, 0, len(noRanks))
		for r := range noRanks {
			noRanksList = append(noRanksList, r)
		}

		if listOrder {
			orders := make([]stringutil.StringCount, 0, len(rankOrder))
			for r, o := range rankOrder {
				orders = append(orders, stringutil.StringCount{Key: r, Count: o})
			}
			sorts.Quicksort(stringutil.ReversedStringCountList{orders})
			preOrder := -1
			for _, order := range orders {
				// fmt.Printf("%d\t%s\n", order.Count, order.Key)
				if order.Count == preOrder {
					fmt.Printf(",%s", order.Key)
				} else {
					if preOrder != -1 {
						fmt.Println()
					}
					fmt.Printf("%s", order.Key)
					preOrder = order.Count
				}
			}
			fmt.Println()
			return
		}

		taxondb := loadTaxonomy(opt, true)

		if opt.Verbose {
			log.Infof("checking defined taxonomic rank order")
		}
		notDefined := make([]string, 0, 10)
		for rank := range taxondb.Ranks {
			if _, ok := rankOrder[rank]; !ok {
				if _, ok := noRanks[rank]; !ok {
					notDefined = append(notDefined, rank)
				}
			}
		}
		if len(notDefined) > 0 {
			checkError(fmt.Errorf("rank order not defined in rank file: %s", strings.Join(notDefined, ", ")))
		}
		if opt.Verbose {
			log.Infof("checking defined taxonomic rank order passed")
		}

		if listRanks {
			orders := make([]stringutil.StringCount, 0, len(taxondb.Ranks))
			var ok bool
			for rank := range taxondb.Ranks {
				if _, ok = rankOrder[rank]; !ok {
					if _, ok := noRanks[rank]; !ok {
						checkError(fmt.Errorf("rank order not defined: %s", rank))
					}
				}
				orders = append(orders, stringutil.StringCount{Key: rank, Count: rankOrder[rank]})
			}
			sorts.Quicksort(stringutil.ReversedStringCountList{orders})
			for _, order := range orders {
				// fmt.Printf("%d\t%s\n", order.Count, order.Key)
				fmt.Printf("%s\n", order.Key)
			}
			return
		}

		tmp := make([]string, 0, len(blackListRanks))
		for _, r := range blackListRanks {
			if r == "" {
				continue
			}
			tmp = append(tmp, strings.ToLower(r))
			blackListRanks = tmp
		}

		if opt.Verbose {
			if discardNoRank {
				log.Infof("ranks without order will be discarded: %s", strings.Join(noRanksList, ", "))
			}
			if saveNorank {
				log.Infof("but some predictable 'no rank' will be saved with -n/--save-predictable-norank)")
			}
			if discardRoot {
				log.Infof("root rank without order will be discarded")
			}
			if len(blackListRanks) > 0 {
				log.Infof("ranks in black list will be discarded: %s", strings.Join(blackListRanks, ", "))
			}
		}

		filter, err := newRankFilter(taxondb, rankOrder, noRanks, lower, higher, equals, blackListRanks, discardNoRank, saveNorank)
		checkError(err)

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

		var writer *unik.Writer

		var infh *bufio.Reader
		var r *os.File
		var reader0 *unik.Reader
		var code uint64
		var taxid uint32
		var k int = -1
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

				reader, err := unik.NewReader(infh)
				checkError(errors.Wrap(err, file))

				hasTaxid = !opt.IgnoreTaxid && reader.HasTaxidInfo()
				if k == -1 {
					reader0 = reader
					k = reader.K
					if !hasTaxid {
						checkError(fmt.Errorf(`taxid information not found: %s`, file))
					}

					mode := reader.Flag
					mode |= unik.UnikIncludeTaxID
					writer, err = unik.NewWriter(outfh, k, mode)
					checkError(errors.Wrap(err, outFile))
					writer.SetMaxTaxid(maxUint32N(reader.GetTaxidBytesLength())) // follow reader
				} else {
					checkCompatibility(reader0, reader, file)
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
						checkError(errors.Wrap(err, file))
					}

					if discardRoot && taxid == rootTaxid {
						continue
					}

					pass, err = filter.isPassed(taxid)
					if err != nil {
						checkError(errors.Wrapf(err, "file: %s, rank: %s", file, rank))
					}
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
	rfilterCmd.Flags().StringP("rank-file", "r", "", `user-defined ordered taxonomic ranks, type "unikmer rfilter --help" for details`)
	rfilterCmd.Flags().BoolP("list-order", "", false, "list defined ranks in order")
	rfilterCmd.Flags().BoolP("list-ranks", "", false, "list ordered ranks in taxonomy database")

	rfilterCmd.Flags().BoolP("discard-noranks", "N", false, `discard ranks without order, type "unikmer filter --help" for details`)
	rfilterCmd.Flags().BoolP("save-predictable-norank", "n", false, `do not discard some special ranks without order when using -L, where rank of the closest higher node is still lower than rank cutoff`)
	rfilterCmd.Flags().StringSliceP("black-list", "B", []string{}, `black list of ranks to discard, e.g., '"no rank", "clade"'`)

	rfilterCmd.Flags().BoolP("discard-root", "R", false, `discard root taxid, defined by --root-taxid`)
	rfilterCmd.Flags().Uint32P("root-taxid", "", 1, `root taxid`)

	rfilterCmd.Flags().StringP("lower-than", "L", "", "output ranks lower than a rank, exclusive with --higher-than")
	rfilterCmd.Flags().StringP("higher-than", "H", "", "output ranks higher than a rank, exclusive with --lower-than")
	rfilterCmd.Flags().StringSliceP("equal-to", "E", []string{}, `output taxIDs with rank equal to some ranks, multiple values can be separated with comma "," (e.g., -E "genus,species"), or give multiple times (e.g., -E genus -E species)`)
}

type rankFilter struct {
	taxondb *taxdump.Taxonomy

	dbRanks   map[string]interface{}
	rankOrder map[string]int

	lower  string
	higher string
	equals []string

	oLower  int
	oHigher int
	oEqual  int

	limitLower  bool
	limitHigher bool
	limitEqual  bool

	noRanks    map[string]interface{}
	blackLists map[string]interface{}
	oEquals    map[int]interface{}

	discardNorank   bool
	saveKnownNoRank bool

	cache map[uint32]bool
}

func newRankFilter(taxondb *taxdump.Taxonomy, rankOrder map[string]int, noRanks map[string]interface{},
	lower string, higher string, equals []string, blackList []string, discardNorank bool, saveKnownNoRank bool) (*rankFilter, error) {

	if lower != "" && higher != "" {
		return nil, fmt.Errorf("higher and lower can't be simultaneous given")
	}

	blackListMap := make(map[string]interface{})
	for _, r := range blackList {
		blackListMap[r] = struct{}{}
	}
	dbRanks := taxondb.Ranks
	f := &rankFilter{
		taxondb:         taxondb,
		dbRanks:         dbRanks,
		rankOrder:       rankOrder,
		lower:           lower,
		higher:          higher,
		equals:          equals,
		noRanks:         noRanks,
		blackLists:      blackListMap,
		discardNorank:   discardNorank,
		saveKnownNoRank: saveKnownNoRank,
		cache:           make(map[uint32]bool, 1024),
	}
	var err error
	if lower != "" {
		f.oLower, err = getRankOrder(dbRanks, rankOrder, lower)
		if err != nil {
			return nil, err
		}
		f.limitLower = true
	}
	if higher != "" {
		f.oHigher, err = getRankOrder(dbRanks, rankOrder, higher)
		if err != nil {
			return nil, err
		}
		f.limitHigher = true
	}
	if len(equals) > 0 {
		f.oEquals = make(map[int]interface{}, len(equals))
		var oe int
		for _, equal := range equals {
			oe, err = getRankOrder(dbRanks, rankOrder, equal)
			if err != nil {
				return nil, err
			}
			f.oEquals[oe] = struct{}{}
		}
		f.limitEqual = true
	}
	return f, nil
}

func getRankOrder(dbRanks map[string]interface{}, rankOrder map[string]int, rank string) (int, error) {
	var ok bool
	if _, ok = rankOrder[rank]; !ok {
		return -1, fmt.Errorf("rank order not defined in rank file: %s", rank)
	}
	if _, ok = dbRanks[rank]; !ok {
		return -1, fmt.Errorf("rank order not found in taxonomy database: %s", rank)
	}

	return rankOrder[rank], nil
}

func (f *rankFilter) isPassed(taxid uint32) (bool, error) {
	rank := f.taxondb.Rank(taxid)
	if rank == "" {
		return false, nil
	}

	rank = strings.ToLower(rank)

	if v, ok := f.cache[taxid]; ok {
		return v, nil
	}

	if _, ok := f.blackLists[rank]; ok {
		f.cache[taxid] = false
		return false, nil
	}

	var isNoRank bool
	_, ok := f.noRanks[rank]
	if ok {
		if f.discardNorank {
			isNoRank = true
			if !f.saveKnownNoRank {
				f.cache[taxid] = false
				return false, nil
			}
		}
	}

	var pass bool

	if isNoRank && f.limitLower && f.saveKnownNoRank {
		nodes := f.taxondb.Nodes
		var _rank string
		var _ok bool
		var _order int

		parent := nodes[taxid]
		for {
			if parent == 1 {
				f.cache[taxid] = false
				return false, nil
			}

			_rank = f.taxondb.Rank(parent)
			_order, _ok = f.rankOrder[_rank]
			if _ok {
				pass = _order <= f.oLower

				f.cache[taxid] = pass
				return pass, nil
			}
			parent = nodes[parent]
		}
	}

	order, _ := f.rankOrder[rank]
	// order, ok := f.rankOrder[rank]
	// if !ok {
	// 	return false, fmt.Errorf("rank order not defined in rank file: %s", rank)
	// }

	if f.limitEqual {
		if _, pass = f.oEquals[order]; pass {
			// pass = true
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

	f.cache[taxid] = pass
	return pass, nil
}

func readRankOrderFromFile(file string) (map[string]int, map[string]interface{}, error) {
	fh, err := os.Open(file)
	if err != nil {
		return nil, nil, fmt.Errorf("read rank order list from '%s': %s", file, err)
	}

	ranks := make([][]string, 0, 128)
	noranks := make(map[string]interface{}, 10)

	scanner := bufio.NewScanner(fh)
	var record, item string
	for scanner.Scan() {
		record = strings.TrimSpace(scanner.Text())
		if record == "" || record[0] == '#' {
			continue
		}

		items := make([]string, 0, 1)

		for _, item = range strings.Split(record, ",") {
			if len(item) == 0 {
				continue
			}
			item = strings.ToLower(strings.TrimSpace(item))

			if item[0] == '!' {
				noranks[item[1:]] = struct{}{}
			} else {
				items = append(items, item)
			}
		}

		if len(items) > 0 {
			ranks = append(ranks, items)
		}
	}
	if err = scanner.Err(); err != nil {
		return nil, nil, fmt.Errorf("read rank order list from '%s': %s", file, err)
	}

	if len(ranks) == 0 {
		return nil, nil, fmt.Errorf("no ranks found in file: %s", file)
	}

	rankOrder := make(map[string]int, len(ranks))
	order := 1
	var ok bool
	var rank string
	for i := len(ranks) - 1; i >= 0; i-- {
		for _, rank = range ranks[i] {
			if _, ok = rankOrder[rank]; ok {
				return nil, nil, fmt.Errorf("duplicated rank: %s", ranks[i])
			}
			rankOrder[rank] = order
		}
		order++
	}
	return rankOrder, noranks, nil
}

func readRankOrder(opt *Options, rankFile string) (map[string]int, map[string]interface{}, error) {
	if rankFile != "" {
		if opt.Verbose {
			log.Infof("read rank order from: %s", rankFile)
		}
		return readRankOrderFromFile(rankFile)
	}

	defaultRankFile := filepath.Join(opt.DataDir, defaultRanksFile)
	existed, err := pathutil.Exists(defaultRankFile)
	if err != nil {
		return nil, nil, fmt.Errorf("check default rank file: %s", defaultRankFile)
	}
	if !existed {
		if opt.Verbose {
			log.Infof("write default rank order to: %s", defaultRankFile)
		}
		err = writeDefaltRankOrderFile(defaultRankFile)
		if err != nil {
			return nil, nil, fmt.Errorf("write default rank file: %s", defaultRankFile)
		}
	}

	if opt.Verbose {
		log.Infof("read rank order from: %s", defaultRankFile)
	}
	return readRankOrderFromFile(defaultRankFile)
}

func writeDefaltRankOrderFile(file string) error {
	return os.WriteFile(file, []byte(defaultRanksText), 0644)
}

const defaultRanksFile = "ranks.txt"
const defaultRanksText = `
# This file defines taxonomic rank order for unikmer/taxonkit.
# 
# Here'are the rules:
#     1. Blank lines or lines starting with "#" are ignored.
#     2. Ranks are in decending order and case ignored.
#     3. Ranks with same order should be in one line separated with comma (",", no space).
#     4. Ranks without order should be assigned a prefix symbol "!" for each rank.
# 
# Deault ranks reference from https://en.wikipedia.org/wiki/Taxonomic_rank ,
# and contains some ranks from NCIB Taxonomy database.
#

!no rank
!clade


life

domain,superkingdom,realm,empire

kingdom
subkingdom
infrakingdom
parvkingdom

superphylum,superdivision
phylum,division
subphylum,subdivision
infraphylum,infradivision
microphylum,microdivision

superclass
class
subclass
infraclass
parvclass

superlegion
legion
sublegion
infralegion

supercohort
cohort
subcohort
infracohort

gigaorder
magnorder,megaorder
grandorder,capaxorder
mirorder,hyperorder
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


superspecies,species group
species subgroup
species

subspecies,forma specialis,pathovar

pathogroup,serogroup
biotype,serotype,genotype

variety,varietas,morph,aberration
subvariety,subvarietas,submorph,subaberration
form,forma
subform,subforma

strain
isolate
`
