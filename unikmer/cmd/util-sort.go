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
	"container/heap"
	"fmt"
	"io"
	"os"
	"path/filepath"

	"github.com/pkg/errors"
	"github.com/shenwei356/bio/taxdump"
	"github.com/shenwei356/unik/v5"
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

	writer, err := unik.NewWriter(outfh, k, mode)
	checkError(err)
	writer.SetMaxTaxid(opt.MaxTaxid)

	var n int64
	var last = ^uint64(0)

	if unique {
		for _, code := range m {
			if code != last {
				writer.WriteCode(code)
				n++
				last = code
			}
		}
	} else if repeated {
		var count int
		var code uint64
		for _, code = range m {
			if code == last {
				count++
				continue
			}

			if count > 0 { // not the first one
				// write all codes once
				writer.WriteCode(last)
				n++
				// write repeated one another time
				if count > 1 {
					writer.WriteCode(last)
					n++
				}
			}

			count = 1
			last = code
		}
		// write all codes once
		writer.WriteCode(last)
		n++
		// write repeated one another time
		if count > 1 {
			writer.WriteCode(last)
			n++
		}
	} else {
		for _, code := range m {

			writer.WriteCode(code)
			n++
		}
	}

	checkError(writer.Flush())
	return n
}

func dumpCodesTaxids2File(mt []CodeTaxid, taxondb *taxdump.Taxonomy, k int, mode uint32, outFile string, opt *Options, unique bool, repeated bool) int64 {
	outfh, gw, w, err := outStream(outFile, opt.Compress, opt.CompressionLevel)
	checkError(err)
	defer func() {
		outfh.Flush()
		if gw != nil {
			gw.Close()
		}
		w.Close()
	}()

	writer, err := unik.NewWriter(outfh, k, mode)
	checkError(err)
	writer.SetMaxTaxid(opt.MaxTaxid)

	var n int64
	var last uint64 = ^uint64(0)
	var lca uint32

	if unique {
		var first bool = true
		for _, codeT := range mt {
			// same k-mer, compute LCA and handle it later
			if codeT.Code == last {
				lca = taxondb.LCA(codeT.Taxid, lca)
				continue
			}

			if first { // just ignore first code, faster than comparing code or slice index, I think
				first = false
			} else { // when meeting new k-mer, output previous one
				writer.WriteCodeWithTaxid(last, lca)
				n++
			}

			last = codeT.Code
			lca = codeT.Taxid
		}
		// do not forget the last one
		writer.WriteCodeWithTaxid(last, lca)
		n++
	} else if repeated {
		var count int
		var codeT CodeTaxid
		for _, codeT = range mt {
			// same k-mer, compute LCA and handle it later
			if codeT.Code == last {
				lca = taxondb.LCA(codeT.Taxid, lca)
				count++
				continue
			}

			if count > 0 { // not the first one
				// write all codes once
				writer.WriteCodeWithTaxid(last, lca)
				n++
				// write repeated one another time
				if count > 1 {
					writer.WriteCodeWithTaxid(last, lca)
					n++
				}
			}

			count = 1
			last = codeT.Code
			lca = codeT.Taxid
		}
		// the last one
		// write all codes once
		writer.WriteCodeWithTaxid(last, lca)
		n++
		// write repeated one another time
		if count > 1 {
			writer.WriteCodeWithTaxid(last, lca)
			n++
		}
	} else {
		writer.Number = uint64(len(mt))
		for _, codeT := range mt {
			writer.WriteCodeWithTaxid(codeT.Code, codeT.Taxid)
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
	idx   int // chunk file index
	code  uint64
	taxid uint32
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

func mergeChunksFile(opt *Options, taxondb *taxdump.Taxonomy, files []string, outFile string, k int, mode uint32, unique bool, repeated bool, finalRound bool) (int64, string) {
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
	hasTaxid := mode&unik.UnikIncludeTaxID > 0
	if hasTaxid && taxondb == nil {
		checkError(fmt.Errorf("taxon information is need when UnikIncludeTaxID is one"))
	}

	writer, err = unik.NewWriter(outfh, k, mode)
	checkError(err)
	writer.SetMaxTaxid(opt.MaxTaxid)

	readers := make(map[int]*unik.Reader, len(files))
	fhs := make([]*os.File, len(files))

	var reader *unik.Reader
	for i, file := range files {
		infh, fh, _, err := inStream(file)
		checkError(errors.Wrap(err, file))
		fhs = append(fhs, fh)

		reader, err := unik.NewReader(infh)
		checkError(errors.Wrap(err, file))
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
		var code uint64
		var taxid uint32
		for i, reader := range readers {
			reader = readers[i]
			n := 0
			for {
				code, taxid, err = reader.ReadCodeWithTaxid()
				if err != nil {
					if err == io.EOF {
						delete(readers, i)
						break
					}
					checkError(fmt.Errorf("faild to fill bufer from file '%s': %s", files[i], err))
				}
				n++
				heap.Push(codes, &codeEntry{idx: i, code: code, taxid: taxid})
				if n >= maxChunkElem {
					break
				}
			}
		}

		return nil
	}

	var e *codeEntry
	var n int64
	var first bool = true
	var last = ^uint64(0)
	var lca uint32
	var code uint64
	var taxid uint32
	var count int

	if hasTaxid {
		if unique {
			for {
				if len(*(codes.entries)) == 0 {
					checkError(fillBuffer())
				}
				if len(*(codes.entries)) == 0 {
					break
				}

				e = heap.Pop(codes).(*codeEntry)
				code = e.code
				taxid = e.taxid

				// -------------------------------------------------

				if code == last {
					lca = taxondb.LCA(taxid, lca)
				} else {
					if first { // just ignore first code, faster than comparing code or slice index, I think
						first = false
					} else { // when meeting new k-mer, output previous one
						writer.WriteCodeWithTaxid(last, lca)
						n++
					}

					last = code
					lca = taxid
				}

				// -------------------------------------------------

				reader = readers[e.idx]
				if reader != nil {
					code, taxid, err = reader.ReadCodeWithTaxid()
					if err != nil {
						if err == io.EOF {
							delete(readers, e.idx)
							continue
						}
						checkError(fmt.Errorf("faild to read from file '%s': %s", files[e.idx], err))
					}
					heap.Push(codes, &codeEntry{idx: e.idx, code: code, taxid: taxid})
				}
			}

			// the last one
			writer.WriteCodeWithTaxid(last, lca)
			n++
		} else if repeated {
			for {
				if len(*(codes.entries)) == 0 {
					checkError(fillBuffer())
				}
				if len(*(codes.entries)) == 0 {
					break
				}

				e = heap.Pop(codes).(*codeEntry)
				code = e.code
				taxid = e.taxid

				// -------------------------------------------------

				// same k-mer, compute LCA and handle it later
				if code == last {
					lca = taxondb.LCA(taxid, lca)
					count++
				} else {
					if count > 0 { // not the first one
						if !finalRound {
							// write all codes once
							writer.WriteCodeWithTaxid(last, lca)
							n++
						}
						// write repeated one another time
						if count > 1 {
							writer.WriteCodeWithTaxid(last, lca)
							n++
						}
					}

					count = 1
					last = code
					lca = taxid
				}

				// -------------------------------------------------

				reader = readers[e.idx]
				if reader != nil {
					code, taxid, err = reader.ReadCodeWithTaxid()
					if err != nil {
						if err == io.EOF {
							delete(readers, e.idx)
							continue
						}
						checkError(fmt.Errorf("faild to read from file '%s': %s", files[e.idx], err))
					}
					heap.Push(codes, &codeEntry{idx: e.idx, code: code, taxid: taxid})
				}
			}

			// the last one
			if count > 0 { // not the first one
				if !finalRound {
					// write all codes once
					writer.WriteCodeWithTaxid(last, lca)
					n++
				}
				// write repeated one another time
				if count > 1 {
					writer.WriteCodeWithTaxid(last, lca)
					n++
				}
			}
		} else {
			for {
				if len(*(codes.entries)) == 0 {
					checkError(fillBuffer())
				}
				if len(*(codes.entries)) == 0 {
					break
				}

				e = heap.Pop(codes).(*codeEntry)
				code = e.code
				taxid = e.taxid

				// -------------------------------------------------

				writer.WriteCodeWithTaxid(code, taxid)
				n++

				// -------------------------------------------------

				reader = readers[e.idx]
				if reader != nil {
					code, taxid, err = reader.ReadCodeWithTaxid()
					if err != nil {
						if err == io.EOF {
							delete(readers, e.idx)
							continue
						}
						checkError(fmt.Errorf("faild to read from file '%s': %s", files[e.idx], err))
					}
					heap.Push(codes, &codeEntry{idx: e.idx, code: code, taxid: taxid})
				}
			}
		}
	} else {
		if unique {
			for {
				if len(*(codes.entries)) == 0 {
					checkError(fillBuffer())
				}
				if len(*(codes.entries)) == 0 {
					break
				}

				e = heap.Pop(codes).(*codeEntry)
				code = e.code
				taxid = e.taxid

				// -------------------------------------------------

				if code != last {
					writer.WriteCode(code)
					n++
					last = code
				}

				// -------------------------------------------------

				reader = readers[e.idx]
				if reader != nil {
					code, taxid, err = reader.ReadCodeWithTaxid()
					if err != nil {
						if err == io.EOF {
							delete(readers, e.idx)
							continue
						}
						checkError(fmt.Errorf("faild to read from file '%s': %s", files[e.idx], err))
					}
					heap.Push(codes, &codeEntry{idx: e.idx, code: code, taxid: taxid})
				}
			}

			// the last one
			if code != last {
				writer.WriteCode(code)
				n++
			}
		} else if repeated {
			for {
				if len(*(codes.entries)) == 0 {
					checkError(fillBuffer())
				}
				if len(*(codes.entries)) == 0 {
					break
				}

				e = heap.Pop(codes).(*codeEntry)
				code = e.code
				taxid = e.taxid

				// -------------------------------------------------

				if code == last {
					count++
				} else {
					if count > 0 { // not the first one
						if !finalRound {
							// write all codes once
							writer.WriteCode(last)
							n++
						}
						// write repeated one another time
						if count > 1 {
							writer.WriteCode(last)
							n++
						}
					}

					count = 1
					last = code
				}

				// -------------------------------------------------

				reader = readers[e.idx]
				if reader != nil {
					code, taxid, err = reader.ReadCodeWithTaxid()
					if err != nil {
						if err == io.EOF {
							delete(readers, e.idx)
							continue
						}
						checkError(fmt.Errorf("faild to read from file '%s': %s", files[e.idx], err))
					}
					heap.Push(codes, &codeEntry{idx: e.idx, code: code, taxid: taxid})
				}
			}

			// the last one
			if code != last {
				if count > 0 { // not the first one
					if !finalRound {
						// write all codes once
						writer.WriteCode(last)
						n++
					}
					// write repeated one another time
					if count > 1 {
						writer.WriteCode(last)
						n++
					}
				}
			}
		} else {
			for {
				if len(*(codes.entries)) == 0 {
					checkError(fillBuffer())
				}
				if len(*(codes.entries)) == 0 {
					break
				}

				e = heap.Pop(codes).(*codeEntry)
				code = e.code
				taxid = e.taxid

				// -------------------------------------------------

				writer.WriteCode(code)
				n++

				// -------------------------------------------------

				reader = readers[e.idx]
				if reader != nil {
					code, taxid, err = reader.ReadCodeWithTaxid()
					if err != nil {
						if err == io.EOF {
							delete(readers, e.idx)
							continue
						}
						checkError(fmt.Errorf("faild to read from file '%s': %s", files[e.idx], err))
					}
					heap.Push(codes, &codeEntry{idx: e.idx, code: code, taxid: taxid})
				}
			}
		}
	}

	checkError(writer.Flush())

	return n, outFile
}
