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
	"compress/flate"
	"fmt"
	"io"
	"sort"
	"strings"

	"github.com/shenwei356/unikmer"
	"github.com/shenwei356/util/pathutil"
	"github.com/spf13/cobra"
)

var mapInitSize = 100000

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
}

func getOptions(cmd *cobra.Command) *Options {
	level := getFlagInt(cmd, "compression-level")
	if level < flate.HuffmanOnly || level > flate.BestCompression {
		checkError(fmt.Errorf("gzip: invalid compression level: %d", level))
	}
	return &Options{
		NumCPUs: getFlagPositiveInt(cmd, "threads"),
		// NumCPUs: 1,
		Verbose:          getFlagBool(cmd, "verbose"),
		Compress:         !getFlagBool(cmd, "no-compress"),
		Compact:          getFlagBool(cmd, "compact"),
		CompressionLevel: level,
	}
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
	dseqs = [][]byte{[]byte{}}
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

func checkFiles(suffix string, files ...string) {
	for _, file := range files {
		if isStdin(file) {
			continue
		}
		ok, err := pathutil.Exists(file)
		if err != nil {
			checkError(fmt.Errorf("fail to read file %s: %s", file, err))
		}
		if !ok {
			checkError(fmt.Errorf("file (linked file) does not exist: %s", file))
		}
		if suffix != "" && !strings.HasSuffix(file, suffix) {
			checkError(fmt.Errorf("input should be stdin or %s file: %s", extDataFile, file))
		}
	}
}

func sortUnikFile(opt Options, unique bool, file string, outFile string) (*unikmer.Header, int, error) {
	// in
	infh, r, _, err := inStream(file)
	if err != nil {
		return nil, 0, err
	}
	defer r.Close()

	var reader *unikmer.Reader
	reader, err = unikmer.NewReader(infh)
	if err != nil {
		return nil, 0, err
	}
	k := reader.K

	// out
	if !isStdout(outFile) {
		outFile += extDataFile
	}
	outfh, gw, w, err := outStream(outFile, opt.Compress, opt.CompressionLevel)
	if err != nil {
		return nil, 0, err
	}
	defer func() {
		outfh.Flush()
		if gw != nil {
			gw.Close()
		}
		w.Close()
	}()

	// mode
	var mode uint32
	if opt.Compact {
		mode |= unikmer.UNIK_COMPACT
	}
	if reader.Flag&unikmer.UNIK_CANONICAL > 0 {
		mode |= unikmer.UNIK_CANONICAL
	}
	mode |= unikmer.UNIK_SORTED

	var writer *unikmer.Writer
	writer, err = unikmer.NewWriter(outfh, k, mode)
	if err != nil {
		return nil, 0, err
	}

	// read
	m := make([]uint64, 0, mapInitSize)
	var kcode unikmer.KmerCode
	for {
		kcode, err = reader.Read()
		if err != nil {
			if err == io.EOF {
				break
			}
			if err != nil {
				return nil, 0, err
			}
		}

		m = append(m, kcode.Code)
	}

	// sort
	sort.Sort(unikmer.CodeSlice(m))

	var n int
	if unique {
		var last = ^uint64(0)
		for _, code := range m {
			if code == last {
				continue
			}
			last = code
			n++
			writer.Write(unikmer.KmerCode{Code: code, K: k})
		}
	} else {
		for _, code := range m {
			writer.Write(unikmer.KmerCode{Code: code, K: k})
		}
		n = len(m)
	}

	return &reader.Header, n, nil
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
