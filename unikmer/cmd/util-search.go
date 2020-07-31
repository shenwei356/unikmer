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
	"os"
	"path/filepath"
	"sync"

	"github.com/edsrzf/mmap-go"
	"github.com/shenwei356/unikmer/index"
	"github.com/shenwei356/util/pathutil"
	"gopkg.in/yaml.v2"
)

const dbInfoFile = "_db.yml"

type UnikIndexDBInfo struct {
	Version   int      `yaml:"version"`
	K         int      `yaml:"k"`
	Canonical bool     `yaml:"canonical"`
	NumHashes int      `yaml:"hashes"`
	FPR       float64  `yaml:"fpr"`
	Kmers     int      `yaml:"totalKmers"`
	Files     []string `yaml:"files"`
	Names     []string `yaml:"names"`
	Sizes     []uint64 `yaml:"kmers"`

	path string
}

func (i UnikIndexDBInfo) String() string {
	return fmt.Sprintf("unikmer index db v%d: k: %d, canonical: %v, #hashes: %d, fpr:%f, #blocks: %d, #%d-mers: %d",
		i.Version, i.K, i.Canonical, i.NumHashes, i.FPR, len(i.Files), i.K, i.Kmers)
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
	return fmt.Sprintf("unikmer index db v%d: #blocks: %d, #%d-mers: %d, #hashes: %d, #signatures: %d",
		db.Info.Version, len(db.Info.Files), db.Header.K, db.Info.Kmers, db.Header.NumHashes, db.Header.NumSigs)
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

func (db *UnikIndexDB) Search(kmers []uint64, threads int, queryCov float64, targetCov float64) map[string][]float64 {
	if threads <= 0 {
		threads = len(db.Indices)
	}
	m := make(map[string][]float64, 8)

	ch := make(chan map[string][]float64, len(db.Indices))
	done := make(chan int)
	go func() {
		for _m := range ch {
			for k, v := range _m {
				m[k] = v
			}
		}
		done <- 1
	}()

	var wg sync.WaitGroup
	tokens := make(chan int, threads)
	for _, idx := range db.Indices {
		wg.Add(1)
		tokens <- 1
		go func(idx *UnikIndex) {
			ch <- idx.Search(kmers, queryCov, targetCov)
			wg.Done()
			<-tokens
		}(idx)
	}

	wg.Wait()
	close(ch)
	<-done

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
	}
	return idx, nil
}

func (idx *UnikIndex) Search(kmers []uint64, queryCov float64, targetCov float64) map[string][]float64 {
	m := make(map[int]int, 8)

	numHashes := int(idx.Header.NumHashes)
	numRowBytes := idx.Header.NumRowBytes
	numSigs := idx.Header.NumSigs
	offset0 := idx.offset0
	fh := idx.fh
	names := idx.Header.Names
	sizes := idx.Header.Sizes
	useMmap := idx.useMmap
	sigs := []byte(idx.sigs)

	data := make([][]uint8, numHashes)
	var row []byte
	var and []byte
	for _, kmer := range kmers {
		for i, loc := range hashLocations(kmer, numHashes, numSigs) {
			row = make([]byte, numRowBytes)

			if useMmap {
				copy(row, sigs[int(offset0+int64(loc*numRowBytes)):int(offset0+int64(loc*numRowBytes))+len(row)])
			} else {
				fh.Seek(offset0+int64(loc*numRowBytes), 0)
				io.ReadFull(fh, row)
			}

			data[i] = row
		}

		// AND
		and = data[0]
		if numHashes > 1 {
			for _, row := range data[1:] {
				for i, b := range row {
					and[i] &= b
				}
			}
		}

		// retrieve info
		for i, b := range and {
			if b == 0 {
				continue
			}
			for j := 0; j < 8; j++ {
				if b&(1<<j) > 0 { // jth bit is 1
					m[i*8+7-j]++
				}
			}
		}
	}

	result := make(map[string][]float64, len(m))
	var t, T float64
	for i, v := range m {
		t = float64(v) / float64(len(kmers))
		if t < queryCov {
			continue
		}
		T = float64(v) / float64(sizes[i])
		if T < targetCov {
			continue
		}

		result[names[i]] = []float64{float64(v), t, T}
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
