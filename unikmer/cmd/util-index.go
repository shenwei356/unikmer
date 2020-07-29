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
	"math"
	"os"
	"sync"

	"github.com/shenwei356/unikmer/index"
)

const extIndex = ".uniki"

func CalcSignatureSize(numElements uint64, numHashes int, falsePositiveRate float64) uint64 {
	ratio := float64(-numHashes) / (math.Log(1 - math.Pow(falsePositiveRate, float64(1/numHashes))))
	return uint64(math.Ceil(float64(numElements) * ratio))
}

type UnikFileInfo struct {
	Path  string
	Name  string
	Kmers int64
}

func (i UnikFileInfo) String() string {
	return fmt.Sprintf("UnikFile{Kmers: %d, Path: %s, Name: %s}", i.Kmers, i.Path, i.Name)
}

type UnikFileInfos []UnikFileInfo

func (l UnikFileInfos) Len() int               { return len(l) }
func (l UnikFileInfos) Less(i int, j int) bool { return l[i].Kmers < l[j].Kmers }
func (l UnikFileInfos) Swap(i int, j int)      { l[i], l[j] = l[j], l[i] }

func MergeUnikIndex(opt *Options, prefix string, files []string, outFile string) error {
	if len(files) == 1 {
		os.Rename(files[0], outFile)
		return nil
	}

	rs := make([]*os.File, len(files))
	defer func() {
		for _, r := range rs {
			r.Close()
		}
	}()

	var header index.Header
	names := make([]string, 0, len(files)*8)

	// retrieve header of the first file, and names in all files
	if opt.Verbose {
		log.Infof("%s checking %d index files", prefix, len(files))
	}
	for i, file := range files {
		infh, r, _, err := inStream(file)
		checkError(err)

		_reader, err := index.NewReader(infh)
		checkError(err)
		if i == 0 {
			header = _reader.Header
		} else if !header.Compatible(_reader.Header) {
			checkError(fmt.Errorf("incompatible index files"))
		}

		names = append(names, _reader.Names...)
		r.Close()
	}
	if opt.Verbose {
		log.Infof("%s # of names: %d", prefix, len(names))
	}

	chs := make([]chan []byte, 0, len(files))

	done := make(chan int)
	go func() {
		outfh, gw, w, err := outStream(outFile, false, -1)
		checkError(err)
		defer func() {
			outfh.Flush()
			if gw != nil {
				gw.Close()
			}
			w.Close()
		}()

		writer, err := index.NewWriter(w, header.K, header.Canonical, header.NumHashes, header.NumSigs, names)
		checkError(err)
		defer func() {
			checkError(writer.Flush())
		}()

		var ok bool
		var closed bool
		first := true
		var nNames int
		var _data []byte
		for {
			var row []byte
			if first {
				row = make([]byte, 0, 8*len(files))
			} else {
				row = make([]byte, 0, nNames)
			}
			for _, ch := range chs {
				_data, ok = <-ch
				if !ok {
					closed = true
				}
				row = append(row, _data...)
			}

			checkError(writer.Write(row))

			if first {
				nNames = len(row)
				first = false
			}
			if closed {
				break
			}
		}

		done <- 1
	}()

	var wg sync.WaitGroup
	for _, file := range files {
		infh, r, _, err := inStream(file)
		checkError(err)
		rs = append(rs, r)

		reader, err := index.NewReader(infh)
		checkError(err)

		ch := make(chan []byte, 8)
		chs = append(chs, ch)

		wg.Add(1)
		go func() {
			for {
				var _data []byte
				_data, err = reader.Read()
				if err != nil {
					if err == io.EOF {
						close(ch)
						r.Close()
						wg.Done()
						break
					}
					checkError(err)
				}
				ch <- _data
			}
		}()
	}

	wg.Wait()
	done <- 1
	return nil
}
