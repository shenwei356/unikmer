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
	"io/ioutil"
	"math"
	"os"
	"path/filepath"
	"sync"

	"github.com/edsrzf/mmap-go"
	"github.com/klauspost/cpuid"
	"github.com/shenwei356/unikmer/index"
	"github.com/shenwei356/util/pathutil"
	"github.com/smallnest/ringbuffer"
	"gopkg.in/yaml.v2"
)

const dbInfoFile = "_db.yml"

type UnikIndexDBInfo struct {
	Version   int      `yaml:"version"`
	K         int      `yaml:"k"`
	Canonical bool     `yaml:"canonical"`
	NumHashes int      `yaml:"hashes"`
	FPR       float64  `yaml:"fpr"`
	BlockSize int      `yaml:"blocksize"`
	Kmers     int      `yaml:"totalKmers"`
	Files     []string `yaml:"files"`
	Names     []string `yaml:"names"`
	Sizes     []uint64 `yaml:"kmers"`

	path string
}

func (i UnikIndexDBInfo) String() string {
	return fmt.Sprintf("unikmer index db v%d: k: %d, canonical: %v, #hashes: %d, fpr:%f, #blocksize: %d, #blocks: %d, #%d-mers: %d",
		i.Version, i.K, i.Canonical, i.NumHashes, i.FPR, i.BlockSize, len(i.Files), i.K, i.Kmers)
}

func NewUnikIndexDBInfo(version int, files []string) UnikIndexDBInfo {
	return UnikIndexDBInfo{Version: version, Files: files}
}

func UnikIndexDBInfoFromFile(file string) (UnikIndexDBInfo, error) {
	info := UnikIndexDBInfo{}

	r, err := os.Open(file)
	if err != nil {
		return info, fmt.Errorf("fail to read unikmer index db info file: %s", file)
	}

	data, err := ioutil.ReadAll(r)
	if err != nil {
		return info, fmt.Errorf("fail to read unikmer index db info file: %s", file)
	}

	err = yaml.Unmarshal(data, &info)
	if err != nil {
		return info, fmt.Errorf("fail to unmarshal unikmer index db info")
	}

	r.Close()

	p, _ := filepath.Abs(file)
	info.path = filepath.Dir(p)
	return info, nil
}

func (i UnikIndexDBInfo) WriteTo(file string) error {
	data, err := yaml.Marshal(i)
	if err != nil {
		return fmt.Errorf("fail to marshal uniker index db info")
	}

	w, err := os.Create(file)
	if err != nil {
		return fmt.Errorf("fail to write unikmer index db info file: %s", file)
	}
	_, err = w.Write(data)
	if err != nil {
		return fmt.Errorf("fail to write unikmer index db info file: %s", file)
	}

	w.Close()
	return nil
}

func (i UnikIndexDBInfo) Check() error {
	for _, file := range i.Files {
		file = filepath.Join(i.path, file)
		ok, err := pathutil.Exists(file)
		if err != nil {
			return fmt.Errorf("error on checking unikmer index file: %s: %s", file, err)
		}
		if !ok {
			return fmt.Errorf("unikmer index file missing: %s", file)
		}
	}
	return nil
}

// ------------------------------------------------------------------

type UnikIndexDB struct {
	path string

	Info   UnikIndexDBInfo
	Header index.Header

	Indices []*UnikIndex
}

func (db *UnikIndexDB) String() string {
	return fmt.Sprintf("unikmer index db v%d: #blocksize: %d, #blocks: %d, #%d-mers: %d, #hashes: %d",
		db.Info.Version, db.Info.BlockSize, len(db.Info.Files), db.Header.K, db.Info.Kmers, db.Header.NumHashes)
}

func NewUnikIndexDB(path string, useMmap bool) (*UnikIndexDB, error) {
	info, err := UnikIndexDBInfoFromFile(filepath.Join(path, dbInfoFile))
	if err != nil {
		return nil, err
	}

	err = info.Check()
	if err != nil {
		return nil, err
	}

	indices := make([]*UnikIndex, 0, len(info.Files))

	// first idx
	idx1, err := NewUnixIndex(filepath.Join(path, info.Files[0]), useMmap)

	if info.Version == int(idx1.Header.Version) &&
		info.K == idx1.Header.K &&
		info.Canonical == idx1.Header.Canonical &&
		info.NumHashes == int(idx1.Header.NumHashes) {
	} else {
		checkError(fmt.Errorf("index files not compatible111"))
	}

	checkError(err)
	indices = append(indices, idx1)

	db := &UnikIndexDB{Info: info, Header: idx1.Header, path: path}

	if len(info.Files) == 1 {
		db.Indices = indices
		return db, nil
	}

	ch := make(chan *UnikIndex, len(info.Files)-1)
	done := make(chan int)
	go func() {
		for idx := range ch {
			indices = append(indices, idx)
		}
		done <- 1
	}()

	var wg sync.WaitGroup
	for _, f := range info.Files[1:] {
		f = filepath.Join(path, f)

		wg.Add(1)
		go func(f string) {
			defer wg.Done()

			idx, err := NewUnixIndex(f, useMmap)
			checkError(err)

			if !idx.Header.Compatible(idx1.Header) {
				checkError(fmt.Errorf("index files not compatible"))
			}

			ch <- idx
		}(f)
	}
	wg.Wait()
	close(ch)
	<-done

	db.Indices = indices
	return db, nil
}

func (db *UnikIndexDB) Close() error {
	var err0 error
	for _, idx := range db.Indices {
		err := idx.Close()
		if err != nil && err0 == nil {
			err0 = err
		}
	}
	return err0
}

func (db *UnikIndexDB) Search(kmers []uint64, threads int, queryCov float64, targetCov float64) map[string][3]float64 {
	if threads <= 0 {
		threads = len(db.Indices)
	}

	m := make(map[string][3]float64, 8)

	numHashes := db.Info.NumHashes
	hashes := make([][]uint64, len(kmers))
	for i, kmer := range kmers {
		hashes[i] = hashValues(kmer, numHashes)
	}

	ch := make([]*map[string][3]float64, len(db.Indices))

	var wg sync.WaitGroup
	// tokens := make(chan int, threads)
	tokens := ringbuffer.New(threads) // ringbufer is faster than channel
	// for i, idx := range db.Indices {
	for i := len(db.Indices) - 1; i >= 0; i-- { // start from bigger files
		idx := db.Indices[i]

		wg.Add(1)
		// tokens <- 1
		tokens.WriteByte(0)
		go func(idx *UnikIndex, ch []*map[string][3]float64, i int) {
			_m := idx.Search(hashes, queryCov, targetCov)
			ch[i] = &_m

			wg.Done()
			// <-tokens
			tokens.ReadByte()
		}(idx, ch, i)
	}
	wg.Wait()

	var k string
	var v [3]float64
	for _, _m := range ch {
		for k, v = range *_m {
			m[k] = v
		}
	}

	return m
}

// ------------------------------------------------------------------

type UnikIndex struct {
	Path   string
	Header index.Header

	fh      *os.File
	reader  *index.Reader
	offset0 int64

	useMmap bool
	sigs    mmap.MMap // mapped sigatures
	sigsB   []byte

	_data  [][]uint8
	buffs  [][]byte
	buffsT [][64]byte // 64 is the cache line size

	pospopcnt func(*[8]int32, []byte)
}

func (idx *UnikIndex) String() string {
	return fmt.Sprintf("%s: %s", idx.Path, idx.Header.String())
}

func NewUnixIndex(file string, useMmap bool) (*UnikIndex, error) {
	fh, err := os.Open(file)
	if err != nil {
		return nil, fmt.Errorf("failed to open unikmer index file: %s", file)
	}

	reader, err := index.NewReader(fh)
	if err != nil {
		return nil, fmt.Errorf("failed to read unikmer index file: %s", file)
	}

	offset, err := fh.Seek(0, 1)
	if err != nil {
		return nil, fmt.Errorf("error on seek unikmer index file: %s", file)
	}

	h := index.Header{}
	h.Version = reader.Version
	h.K = reader.K
	h.Canonical = reader.Canonical
	h.NumHashes = reader.NumHashes
	h.Names = reader.Names
	h.Sizes = reader.Sizes
	h.NumRowBytes = reader.NumRowBytes
	h.NumSigs = reader.NumSigs
	idx := &UnikIndex{Path: file, Header: h, fh: fh, reader: reader, offset0: offset}
	idx.useMmap = useMmap
	if useMmap {
		idx.sigs, err = mmap.Map(fh, mmap.RDONLY, 0)
		if err != nil {
			return nil, err
		}
		idx.sigsB = []byte(idx.sigs)
	}
	idx._data = make([][]uint8, reader.NumHashes)
	for i := 0; i < int(reader.NumHashes); i++ {
		idx._data[i] = make([]byte, reader.NumRowBytes)
	}

	// byte matrix for counting
	buffs := make([][]byte, 64)
	for i := 0; i < 64; i++ {
		buffs[i] = make([]byte, reader.NumRowBytes)
	}

	// transpose of buffs
	buffsT := make([][64]byte, reader.NumRowBytes)
	for i := 0; i < reader.NumRowBytes; i++ {
		buffsT[i] = [64]byte{}
	}
	idx.buffs = buffs
	idx.buffsT = buffsT

	if AVX2Available {
		idx.pospopcnt = Pospopcnt
	} else {
		idx.pospopcnt = PospopcntGo
	}
	return idx, nil
}

func (idx *UnikIndex) Search(hashes [][]uint64, queryCov float64, targetCov float64) map[string][3]float64 {
	numNames := len(idx.Header.Names)

	numHashes := int(idx.Header.NumHashes)
	numRowBytes := idx.Header.NumRowBytes
	numSigs := idx.Header.NumSigs
	numSigsInt := uint64(numSigs)
	offset0 := idx.offset0
	fh := idx.fh
	names := idx.Header.Names
	sizes := idx.Header.Sizes
	sigs := idx.sigsB
	useMmap := idx.useMmap
	data := idx._data
	pospopcnt := idx.pospopcnt

	counts := make([][8]int32, numRowBytes)

	var offset int
	var offset2 int64
	var loc int
	var i, j int
	var hs []uint64
	var row []byte
	var b byte
	var h uint64

	buffs := idx.buffs
	buffsT := idx.buffsT
	bufIdx := 0
	var buf *[64]byte

	for _, hs = range hashes {
		if useMmap {
			for i, h = range hs {
				loc = int(h % numSigsInt)
				offset = int(offset0 + int64(loc*numRowBytes))

				data[i] = sigs[offset : offset+numRowBytes]
			}
		} else {
			for i, h = range hs {
				loc = int(h % numSigsInt)
				offset2 = offset0 + int64(loc*numRowBytes)

				fh.Seek(offset2, 0)
				io.ReadFull(fh, data[i])
			}
		}

		// AND
		var and []byte

		and = data[0]
		if numHashes > 1 {
			for _, row = range data[1:] {
				for i, b = range row {
					and[i] &= b
				}
			}
		}

		// add to buffer for counting
		buffs[bufIdx] = and
		bufIdx++

		if bufIdx == 64 {
			// transpose
			for i = 0; i < numRowBytes; i++ { // every column in matrix
				buf = &buffsT[i]
				for j = 0; j < 64; j++ {
					(*buf)[j] = buffs[j][i]
				}
			}

			// count
			for i = 0; i < numRowBytes; i++ { // every row in transposed matrix
				pospopcnt(&counts[i], buffsT[i][:])
			}

			bufIdx = 0
		}
	}
	// left data in buffer
	if bufIdx > 0 {
		// transpose
		for i = 0; i < numRowBytes; i++ { // every column in matrix
			buf = &buffsT[i]
			for j = 0; j < 64; j++ {
				(*buf)[j] = buffs[j][i]
			}
		}

		// count
		for i = 0; i < numRowBytes; i++ { // every row in transposed matrix
			pospopcnt(&counts[i], buffsT[i][0:bufIdx])
		}
	}

	var _counts [8]int32
	var count int32
	var ix8 int
	var k int

	m := make([]int, numNames)

	iLast := numRowBytes - 1
	nTargets := 0
	for i, _counts = range counts {
		ix8 = i << 3
		for j, count = range _counts {
			k = ix8 + j
			if i == iLast && k == numNames {
				break
			}
			if count > 0 {
				nTargets++
			}
			m[k] = int(count)
		}
	}

	result := make(map[string][3]float64, nTargets)
	var t, T float64
	for i, v := range m {
		if v == 0 {
			continue
		}

		t = float64(v) / float64(len(hashes))
		if t < queryCov {
			continue
		}
		T = float64(v) / float64(sizes[i])
		if T < targetCov {
			continue
		}

		result[names[i]] = [3]float64{float64(v), t, T}
	}

	return result
}

func (idx *UnikIndex) Close() error {
	if idx.useMmap {
		err := idx.sigs.Unmap()
		if err != nil {
			return err
		}
	}
	return idx.fh.Close()
}

func maxFPR(p float64, k float64, l int) float64 {
	return math.Exp(-float64(l) * (k - p) * (k - p) / 2 / (1 - p))
}

var AVX2Available = cpuid.CPU.AVX2()

// for CPUs not supporting AVX2
func PospopcntGo(counts *[8]int32, buf []byte) {
	for i := 0; i < len(buf); i++ {
		for j := 0; j < 8; j++ {
			(*counts)[7-j] += int32(buf[i]) >> j & 1
		}
	}
}
