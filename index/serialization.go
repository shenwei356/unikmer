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

package index

import (
	"encoding/binary"
	"errors"
	"fmt"
	"io"
	"strings"
)

// Version is the version of index format
const Version uint8 = 2

// Magic number of index file.
var Magic = [8]byte{'.', 'u', 'n', 'i', 'k', 'i', 'd', 'x'}

// ErrInvalidIndexFileFormat means invalid index format.
var ErrInvalidIndexFileFormat = errors.New("unikmer/index: invalid index format")

// ErrUnfishedWrite means writing not finished
var ErrUnfishedWrite = errors.New("unikmer/index: index not fished writing")

// ErrTruncateIndexFile means the file is truncated
var ErrTruncateIndexFile = errors.New("unikmer/index: truncated index file")

// ErrWrongWriteDataSize means the size of data to write is invalid
var ErrWrongWriteDataSize = errors.New("unikmer/index: write data with wrong size")

// ErrVersionMismatch means version mismatch between files and program
var ErrVersionMismatch = errors.New("unikmer/index: version mismatch")

// ErrNameAndSizeMismatch means size of names and sizes are not equal.
var ErrNameAndSizeMismatch = errors.New("unikmer/index: size of names and sizes unequal")

// ErrNameAndIndexMismatch means size of names and sizes are not equal.
var ErrNameAndIndexMismatch = errors.New("unikmer/index: size of names and indices unequal")

var be = binary.BigEndian

// Header contains metadata
type Header struct {
	Version   uint8 // uint8
	K         int   // uint8
	Canonical bool  // uint8
	NumHashes uint8 // uint8
	NumSigs   uint64
	Names     []string
	Indices   []uint32
	Sizes     []uint64

	NumRowBytes int // length of bytes for storing one row of signiture for n names
}

func (h Header) String() string {
	return fmt.Sprintf("unikmer index file v%d: k: %d, canonical: %v, #hashes: %d, #signatures: %d, names: %s",
		h.Version, h.K, h.Canonical, h.NumHashes, h.NumSigs, strings.Join(h.Names, ", "))
}

// Compatible checks compatibility
func (h Header) Compatible(b Header) bool {
	if h.Version == b.Version &&
		h.K == b.K &&
		h.Canonical == b.Canonical &&
		h.NumHashes == b.NumHashes {

		return true
	}
	return false
}

// Reader is for reading KmerCode.
type Reader struct {
	Header
	r io.Reader

	count uint64
}

// NewReader returns a Reader.
func NewReader(r io.Reader) (reader *Reader, err error) {
	reader = &Reader{r: r}
	err = reader.readHeader()
	if err != nil {
		return nil, err
	}

	reader.NumRowBytes = int((len(reader.Names) + 7) / 8)
	return reader, nil
}

func (reader *Reader) readHeader() (err error) {
	// check Magic number
	var m [8]byte
	r := reader.r
	err = binary.Read(r, be, &m)
	if err != nil {
		return err
	}
	same := true
	for i := 0; i < 8; i++ {
		if Magic[i] != m[i] {
			same = false
			break
		}
	}
	if !same {
		return ErrInvalidIndexFileFormat
	}

	// 4 bytes meta info
	var meta [4]uint8
	err = binary.Read(r, be, &meta)
	if err != nil {
		return err
	}
	// check compatibility
	if Version != meta[0] {
		return ErrVersionMismatch
	}
	reader.Version = meta[0]
	reader.K = int(meta[1])
	if meta[2] > 0 {
		reader.Canonical = true
	}
	reader.NumHashes = meta[3]

	// 8 bytes signature size
	err = binary.Read(r, be, &reader.NumSigs)
	if err != nil {
		return err
	}

	// 4 bytes length of Names
	var n uint32
	err = binary.Read(r, be, &n)
	if err != nil {
		return err
	}

	// Names
	namesData := make([]byte, n)
	err = binary.Read(r, be, &namesData)
	if err != nil {
		return err
	}
	names := strings.Split(string(namesData), "\n")
	names = names[0 : len(names)-1]
	reader.Names = names

	// Indices
	indicesData := make([]uint32, len(names))
	err = binary.Read(r, be, &indicesData)
	if err != nil {
		return err
	}
	reader.Indices = indicesData

	// Sizes
	sizesData := make([]uint64, len(names))
	err = binary.Read(r, be, &sizesData)
	if err != nil {
		return err
	}
	reader.Sizes = sizesData

	return nil
}

// Read reads one code.
func (reader *Reader) Read() ([]byte, error) {
	data := make([]byte, reader.NumRowBytes)
	nReaded, err := io.ReadFull(reader.r, data)
	if err != nil {
		if err == io.EOF {
			if reader.count != reader.NumSigs {
				return nil, ErrTruncateIndexFile
			}
		}
		return nil, err
	}
	if nReaded < reader.NumRowBytes {
		return nil, ErrTruncateIndexFile
	}
	reader.count++
	return data, nil
}

// Writer writes KmerCode.
type Writer struct {
	Header
	w           io.Writer
	wroteHeader bool

	count uint64
}

// NewWriter creates a Writer.
func NewWriter(w io.Writer, k int, canonical bool, numHashes uint8, numSigs uint64, names []string, indices []uint32, sizes []uint64) (*Writer, error) {
	if len(names) != len(sizes) {
		return nil, ErrNameAndSizeMismatch
	}
	if len(names) != len(indices) {
		return nil, ErrNameAndIndexMismatch
	}

	writer := &Writer{
		Header: Header{
			Version:   Version,
			K:         k,
			Canonical: canonical,
			NumHashes: numHashes,
			NumSigs:   numSigs,
			Names:     names,
			Indices:   indices,
			Sizes:     sizes,
		},
		w: w,
	}
	writer.NumRowBytes = int((len(names) + 7) / 8)

	return writer, nil
}

// WriteHeader writes file header
func (writer *Writer) WriteHeader() (err error) {
	if writer.wroteHeader {
		return nil
	}
	w := writer.w

	// 8 bytes magic number
	err = binary.Write(w, be, Magic)
	if err != nil {
		return err
	}

	// 4 bytes meta info
	var canonical uint8
	if writer.Canonical {
		canonical = 1
	}
	err = binary.Write(w, be, [4]uint8{writer.Version, uint8(writer.K), canonical, writer.NumHashes})
	if err != nil {
		return err
	}

	// 8 bytes signature size
	err = binary.Write(w, be, writer.NumSigs)
	if err != nil {
		return err
	}

	// 4 bytes length of Names
	var n int
	for _, name := range writer.Names {
		n += len(name) + 1
	}

	err = binary.Write(w, be, uint32(n))
	if err != nil {
		return err
	}

	// Names
	for _, name := range writer.Names {
		err = binary.Write(w, be, []byte(name+"\n"))
		if err != nil {
			return err
		}
	}

	// Indices
	err = binary.Write(w, be, writer.Indices)
	if err != nil {
		return err
	}

	// Sizes
	err = binary.Write(w, be, writer.Sizes)
	if err != nil {
		return err
	}

	writer.wroteHeader = true
	return nil
}

// Write writes some thing
func (writer *Writer) Write(data []byte) (err error) {
	if len(data) != writer.NumRowBytes {
		return ErrWrongWriteDataSize
	}

	// lazily write header
	if !writer.wroteHeader {
		err = writer.WriteHeader()
		if err != nil {
			return err
		}
		writer.wroteHeader = true
	}

	_, err = writer.w.Write(data)
	if err != nil {
		return err
	}

	writer.count++
	return nil
}

// WriteBatch writes a batch of data
func (writer *Writer) WriteBatch(data []byte, n int) (err error) {
	// lazily write header
	if !writer.wroteHeader {
		err = writer.WriteHeader()
		if err != nil {
			return err
		}
		writer.wroteHeader = true
	}

	_, err = writer.w.Write(data)
	if err != nil {
		return err
	}

	writer.count += uint64(n)
	return nil
}

// Flush check completeness
func (writer *Writer) Flush() (err error) {
	if !writer.wroteHeader {
		writer.WriteHeader()
	}
	if writer.count != writer.NumSigs {
		return ErrUnfishedWrite
	}
	return nil
}
