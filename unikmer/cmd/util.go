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
	"compress/flate"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"strconv"
	"strings"

	"github.com/pkg/errors"
	"github.com/shenwei356/bio/seqio/fastx"
	"github.com/shenwei356/unikmer"
	"github.com/shenwei356/util/pathutil"
	"github.com/spf13/cobra"
	"github.com/twotwotwo/sorts"
)

var mapInitSize = 1 << 20 // 1M

const (
	flagContinue = iota
	flagBreak
	flagReturn
)

// Options contains the global flags
type Options struct {
	NumCPUs          int
	Verbose          bool
	Compress         bool
	Compact          bool
	CompressionLevel int
	MaxTaxid         uint32
	IgnoreTaxid      bool
	DataDir          string
	NodesFile        string
	CacheLCA         bool

	NoCheckFile bool
}

func getOptions(cmd *cobra.Command) *Options {
	level := getFlagInt(cmd, "compression-level")
	if level < flate.HuffmanOnly || level > flate.BestCompression {
		checkError(fmt.Errorf("gzip: invalid compression level: %d", level))
	}

	var val, dataDir string
	if val = os.Getenv("UNIKMER_DB"); val != "" {
		dataDir = val
	} else {
		dataDir = getFlagString(cmd, "data-dir")
	}

	threads := getFlagPositiveInt(cmd, "threads")
	if threads >= 1000 {
		checkError(fmt.Errorf("are your seriously? %d threads? It will exhaust your RAM", threads))
	}

	sorts.MaxProcs = threads

	return &Options{
		NumCPUs:          threads,
		Verbose:          getFlagBool(cmd, "verbose"),
		Compress:         !getFlagBool(cmd, "no-compress"),
		Compact:          getFlagBool(cmd, "compact"),
		CompressionLevel: level,

		MaxTaxid:    getFlagUint32(cmd, "max-taxid"),
		IgnoreTaxid: getFlagBool(cmd, "ignore-taxid"),

		DataDir:  dataDir,
		CacheLCA: true, // getFlagBool(cmd, "cache-lca"),

		NoCheckFile: getFlagBool(cmd, "nocheck-file"),
	}
}

func checkDataDir(opt *Options) {
	existed, err := pathutil.DirExists(opt.DataDir)
	checkError(errors.Wrap(err, opt.DataDir))
	if !existed {
		log.Errorf(`data directory not created. please download and decompress ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz, and copy "nodes.dmp" to %s`, opt.DataDir)
	}
}

func loadTaxonomy(opt *Options, withRank bool) *unikmer.Taxonomy {
	checkDataDir(opt)

	if opt.Verbose {
		log.Infof("loading Taxonomy from: %s", opt.DataDir)
	}
	var t *unikmer.Taxonomy
	var err error
	if withRank {
		t, err = unikmer.NewTaxonomyWithRankFromNCBI(filepath.Join(opt.DataDir, "nodes.dmp"))
	} else {
		t, err = unikmer.NewTaxonomyFromNCBI(filepath.Join(opt.DataDir, "nodes.dmp"))
	}
	if err != nil {
		checkError(fmt.Errorf("err on loading Taxonomy nodes: %s", err))
	}
	if opt.Verbose {
		if withRank {
			log.Infof("%d nodes in %d ranks loaded", len(t.Nodes), len(t.Ranks))
		} else {
			log.Infof("%d nodes loaded", len(t.Nodes))
		}
	}

	var existed bool

	// err = t.LoadDeletedNodesFromNCBI(filepath.Join(opt.DataDir, "delnodes.dmp"))
	// if err != nil {
	// 	checkError(fmt.Errorf("err on loading Taxonomy nodes: %s", err))
	// }

	existed, err = pathutil.Exists(filepath.Join(opt.DataDir, "merged.dmp"))
	if err != nil {
		checkError(fmt.Errorf("err on checking file merged.dmp: %s", err))
	}
	if existed {
		err = t.LoadMergedNodesFromNCBI(filepath.Join(opt.DataDir, "merged.dmp"))
		if err != nil {
			checkError(fmt.Errorf("err on loading Taxonomy merged nodes: %s", err))
		}
	}

	if opt.Verbose {
		log.Infof("%d merged nodes loaded", len(t.MergeNodes))
	}

	if opt.CacheLCA {
		t.CacheLCA()
	}

	opt.MaxTaxid = t.MaxTaxid()
	return t
}

var degenerateBaseMapNucl = map[byte]string{
	'A': "A",
	'T': "T",
	'U': "U",
	'C': "C",
	'G': "G",
	'R': "AG",
	'Y': "CT",
	'M': "AC",
	'K': "GT",
	'S': "CG",
	'W': "AT",
	'H': "ACT",
	'B': "CGT",
	'V': "ACG",
	'D': "AGT",
	'N': "ACGT",
	'a': "a",
	't': "t",
	'u': "u",
	'c': "c",
	'g': "g",
	'r': "ag",
	'y': "ct",
	'm': "ac",
	'k': "gt",
	's': "cg",
	'w': "at",
	'h': "act",
	'b': "cgt",
	'v': "acg",
	'd': "agt",
	'n': "acgt",
}

func extendDegenerateSeq(s []byte) (dseqs [][]byte, err error) {
	dseqs = [][]byte{{}}
	var i, j, k int
	var ok bool
	var dbases string
	var dbase byte
	for _, base := range s {
		if dbases, ok = degenerateBaseMapNucl[base]; ok {
			if len(dbases) == 1 {
				dbase = dbases[0]
				for i = 0; i < len(dseqs); i++ {
					dseqs[i] = append(dseqs[i], dbase)
				}
			} else {
				// 2nd
				more := make([][]byte, len(dseqs)*(len(dbases)-1))
				k = 0
				for i = 1; i < len(dbases); i++ {
					for j = 0; j < len(dseqs); j++ {
						more[k] = []byte(string(append(dseqs[j], dbases[i])))
						k++
					}
				}

				// 1th
				for i = 0; i < len(dseqs); i++ {
					dseqs[i] = append(dseqs[i], dbases[0])
				}

				dseqs = append(dseqs, more...)
			}

		} else {
			return dseqs, unikmer.ErrIllegalBase
		}
	}
	return dseqs, nil
}

func checkFileSuffix(opt *Options, suffix string, files ...string) {
	if opt.NoCheckFile {
		return
	}

	for _, file := range files {
		if isStdin(file) {
			continue
		}

		if suffix != "" && !strings.HasSuffix(file, suffix) {
			checkError(fmt.Errorf("input should be stdin or %s file: %s", suffix, file))
		}
	}
}

func uniqInts(data []int) []int {
	if len(data) == 0 || len(data) == 1 {
		return data
	}
	m := make(map[int]struct{}, len(data))
	for _, d := range data {
		m[d] = struct{}{}
	}
	data2 := make([]int, len(m))
	i := 0
	for k := range m {
		data2[i] = k
		i++
	}
	return data2
}

func minInt(a int, vals ...int) int {
	min := a
	for _, v := range vals {
		if v < min {
			min = v
		}
	}
	return min
}

// ParseByteSize parses byte size from string.
func ParseByteSize(val string) (int, error) {
	val = strings.Trim(val, " \t\r\n")
	if val == "" {
		return 0, nil
	}
	var u int64
	var noUnit bool
	switch val[len(val)-1] {
	case 'B', 'b':
		u = 1
	case 'K', 'k':
		u = 1 << 10
	case 'M', 'm':
		u = 1 << 20
	case 'G', 'g':
		u = 1 << 30
	default:
		noUnit = true
		u = 1
	}
	var size float64
	var err error
	if noUnit {
		size, err = strconv.ParseFloat(val, 10)
		if err != nil {
			return 0, fmt.Errorf("invalid byte size: %s", val)
		}
		if size < 0 {
			size = 0
		}
		return int(size), nil
	}

	if len(val) == 1 { // no value
		return 0, nil
	}

	size, err = strconv.ParseFloat(strings.Trim(val[0:len(val)-1], " \t\r\n"), 10)
	if err != nil {
		return 0, fmt.Errorf("invalid byte size: %s", val)
	}
	if size < 0 {
		size = 0
	}
	return int(size * float64(u)), nil
}

var maxUint32 = uint64(^uint32(0))

func maxUint32N(n int) uint32 {
	return (1 << (n << 3)) - 1
}

func loadHash2Loc(files []string, k int) ([][]byte, map[uint64][2]int, error) {
	var hash2loc map[uint64][2]int // hash -> [seq idx, seq loc]
	var sequences [][]byte

	sequences = make([][]byte, 0, 8)
	hash2loc = make(map[uint64][2]int, mapInitSize)

	var err error
	var fastxReader *fastx.Reader
	var record *fastx.Record
	var iter *unikmer.Iterator
	var code uint64
	var ok bool
	var seqIdx int

	for _, file := range files {
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
			iter, err = unikmer.NewHashIterator(record.Seq, k, true, true)
			if err != nil {
				checkError(errors.Wrapf(err, "seq: %s", record.Name))
			}

			sequences = append(sequences, record.Seq.Clone().Seq)

			for {
				code, ok = iter.NextHash()
				if !ok {
					break
				}

				if _, ok = hash2loc[code]; !ok {
					hash2loc[code] = [2]int{seqIdx, iter.Index()}
				}
			}
			seqIdx++
		}

	}
	return sequences, hash2loc, nil
}
