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
	"container/heap"
	"fmt"
	"io"
	"os"
	"path/filepath"

	"github.com/shenwei356/unikmer"
)

func dumpCodes2File(m []uint64, k int, mode uint32, outFile string, opt *Options, unique bool, repeated bool) int64 {
	outfh, gw, w, err := outStream(outFile, opt.Compress, opt.CompressionLevel)
	checkError(err)
	defer func() {
		outfh.Flush()
		if gw != nil {
			gw.Close()
		}
		w.Close()
	}()

	writer, err := unikmer.NewWriter(outfh, k, mode)
	checkError(err)

	var n int64
	var last = ^uint64(0)
	var count int
	// log.Warningf("%d", m)
	for _, code := range m {
		if unique {
			if code != last {
				writer.WriteCode(code)
				n++
				last = code
			}
		} else if repeated {
			// log.Warningf("last: %d, code: %d, count: %d", last, code, count)
			if code == last {
				if count == 1 { // write once
					writer.WriteCode(code)
					n++
					count++
				}
			} else {
				writer.WriteCode(code)
				n++

				last = code
				count = 1
			}
		} else {
			writer.WriteCode(code)
			n++
		}
	}

	checkError(writer.Flush())
	return n
}

func chunkFileName(outDir string, i int) string {
	return filepath.Join(outDir, fmt.Sprintf("chunk_%03d", i)) + extDataFile
}

type codeEntry struct {
	idx  int // chunk file index
	code uint64
}

type codeEntryHeap struct {
	entries *[]*codeEntry
}

func (h codeEntryHeap) Len() int { return len(*(h.entries)) }

func (h codeEntryHeap) Less(i, j int) bool {
	return (*(h.entries))[i].code < (*(h.entries))[j].code
}

func (h codeEntryHeap) Swap(i, j int) {
	(*(h.entries))[i], (*(h.entries))[j] = (*(h.entries))[j], (*(h.entries))[i]
}

func (h codeEntryHeap) Push(x interface{}) {
	*(h.entries) = append(*(h.entries), x.(*codeEntry))
}

func (h codeEntryHeap) Pop() interface{} {
	n := len(*(h.entries))
	x := (*(h.entries))[n-1]
	*(h.entries) = (*(h.entries))[:n-1]
	return x
}

func mergeChunksFile(opt *Options, files []string, outFile string, k int, mode uint32, unique bool, repeated bool, finalRound bool) (int64, string) {
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
	var kcode unikmer.KmerCode

	writer, err = unikmer.NewWriter(outfh, k, mode)
	checkError(err)

	readers := make(map[int]*unikmer.Reader, len(files))
	fhs := make([]*os.File, len(files))

	var reader *unikmer.Reader
	for i, file := range files {
		infh, fh, _, err := inStream(file)
		checkError(err)
		fhs = append(fhs, fh)

		reader, err := unikmer.NewReader(infh)
		checkError(err)
		readers[i] = reader
	}
	defer func() {
		for _, fh := range fhs {
			fh.Close()
		}
	}()

	maxChunkElem := 2

	entries := make([]*codeEntry, 0, len(files)*maxChunkElem)
	codes := codeEntryHeap{entries: &entries}

	fillBuffer := func() error {
		var err error
		var kcode unikmer.KmerCode
		for i, reader := range readers {
			reader = readers[i]
			n := 0
			for {
				kcode, err = reader.Read()
				if err != nil {
					if err == io.EOF {
						delete(readers, i)
						break
					}
					checkError(fmt.Errorf("faild to fill bufer from file '%s': %s", files[i], err))
				}
				n++
				heap.Push(codes, &codeEntry{idx: i, code: kcode.Code})
				if n >= maxChunkElem {
					break
				}
			}
		}

		return nil
	}

	var e *codeEntry
	var n int64
	var last = ^uint64(0)
	var code uint64
	var count int

	for {
		if len(*(codes.entries)) == 0 {
			checkError(fillBuffer())
		}
		if len(*(codes.entries)) == 0 {
			break
		}

		e = heap.Pop(codes).(*codeEntry)
		code = e.code

		if unique {
			if code != last {
				writer.WriteCode(code)
				n++
				last = code
			}
		} else if repeated {
			if code == last {
				if count == 1 { // write another copy
					writer.WriteCode(code)
					n++
					count++
				}
			} else {
				if !finalRound {
					writer.WriteCode(code)
					n++
				}

				last = code
				count = 1
			}
		} else {
			writer.WriteCode(code)
			n++
		}

		reader = readers[e.idx]
		if reader != nil {
			kcode, err = reader.Read()
			if err != nil {
				if err == io.EOF {
					delete(readers, e.idx)
					continue
				}
				checkError(fmt.Errorf("faild to read from file '%s': %s", files[e.idx], err))
			}
			heap.Push(codes, &codeEntry{e.idx, kcode.Code})
		}
	}

	checkError(writer.Flush())

	return n, outFile
}
