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

package unikmer

import (
	"encoding/binary"
	"errors"
	"fmt"
	"io"
)

// MainVersion is the main version number.
const MainVersion uint8 = 1

// MinorVersion is the minor version number.
const MinorVersion uint8 = 0

// Magic number of binary file.
var Magic = [8]byte{'.', 'u', 'n', 'i', 'k', 'm', 'e', 'r'}

// ErrInvalidFileFormat means invalid file format.
var ErrInvalidFileFormat = errors.New("unikmer: invalid binary format")

// ErrBrokenFile means the file is not complete.
// var ErrBrokenFile = errors.New("unikmer: broken file")

// ErrKMismatch means K size mismatch.
var ErrKMismatch = errors.New("unikmer: K mismatch")

var be = binary.BigEndian

// Header contains metadata
type Header struct {
	MainVersion  uint8
	MinorVersion uint8
	K            int
	Flag         uint32
}

const (
	UNIK_COMPACT   = 1 << iota
	UNIK_CANONICAL = 1 << iota
)

func (h Header) String() string {
	return fmt.Sprintf("unikmer binary kmer data file v%d.%d with K=%d and Flag=%d",
		h.MainVersion, h.MinorVersion, h.K, h.Flag)
}

// Reader is for reading KmerCode.
type Reader struct {
	Header
	r    io.Reader
	err  error
	code uint64
	size uint64

	compact bool // saving KmerCode in variable-length byte array.
	buf     []byte
	bufsize int
}

// NewReader returns a Reader.
func NewReader(r io.Reader) (*Reader, error) {
	reader := &Reader{r: r}
	reader.err = reader.readHeader()
	if reader.err != nil {
		return nil, reader.err
	}
	return reader, nil
}

func (reader *Reader) readHeader() error {
	// check Magic number
	var m [8]byte
	reader.err = binary.Read(reader.r, be, &m)
	if reader.err != nil {
		return reader.err
	}
	same := true
	for i := 0; i < 8; i++ {
		if Magic[i] != m[i] {
			same = false
			break
		}
	}
	if !same {
		return ErrInvalidFileFormat
	}

	// read metadata
	var meta [4]uint8
	reader.err = binary.Read(reader.r, be, &meta)
	if reader.err != nil {
		return reader.err
	}
	// check compatibility？
	if (meta[0] == 0 && meta[1] == 0) ||
		MainVersion != meta[0] {
		return fmt.Errorf(".unik format compatibility error, please recreate with newest version")
	}
	reader.MainVersion = meta[0]
	reader.MinorVersion = meta[1]

	reader.K = int(meta[2])

	reader.buf = make([]byte, 8)
	reader.bufsize = int((reader.K + 3) / 4)

	reader.err = binary.Read(reader.r, be, &reader.Flag)
	if reader.err != nil {
		return reader.err
	}
	if reader.Flag&UNIK_COMPACT > 0 {
		reader.compact = true
	}

	return nil
}

// Read reads one KmerCode.
func (reader *Reader) Read() (KmerCode, error) {
	if reader.compact {
		reader.err = binary.Read(reader.r, be, reader.buf[8-reader.bufsize:])
		reader.code = be.Uint64(reader.buf)
	} else {
		reader.err = binary.Read(reader.r, be, &reader.code)
	}
	if reader.err != nil {
		return KmerCode{}, reader.err
	}

	reader.size++
	return KmerCode{Code: reader.code, K: reader.Header.K}, nil
}

// Writer writes KmerCode.
//
//     # MainVersion: 1
//     offset      bytes   name                          type
//     0           8       magic number(".unikmer")      `[8]byte`
//     64          1       MainVersion                   `uint8`
//     72          1       MinorVersion                  `uint8`
//     80          1       K                             `uint8`
//     88          1       reserved                      `uint8`
//     96          4       Flag                          `uint32`
//     128+i×64    8       (kmers)<sup>a</sup>           `uint64`
//     128+i×8×n   8×n     (compact kmers)<sup>b</sup>   `[n]byte`
//
// - <sup>a</sup> One Kmer is encoded as `uint64` and serialized in 8 Bytes by default.
//
// - <sup>b</sup> In compact mode, `x = int((k + 3) / 4)`.
//
type Writer struct {
	Header
	w           io.Writer
	kcode       KmerCode
	wroteHeader bool
	err         error
	size        int64

	compact bool // saving KmerCode in variable-length byte array.
	buf     []byte
	bufsize int
}

// NewWriter creates a Writer.
func NewWriter(w io.Writer, k int, flag uint32) (*Writer, error) {
	if k == 0 || k > 32 {
		return nil, ErrKOverflow
	}

	writer := &Writer{
		Header:  Header{MainVersion: MainVersion, MinorVersion: MinorVersion, K: k, Flag: flag},
		w:       w,
		buf:     make([]byte, 8),
		bufsize: int(k+3) / 4,
	}
	if writer.Flag&UNIK_COMPACT > 0 {
		writer.compact = true
	}
	return writer, nil
}

func (writer *Writer) writeHeader() error {
	// write magic number
	writer.err = binary.Write(writer.w, be, Magic)
	if writer.err != nil {
		return writer.err
	}

	writer.err = binary.Write(writer.w, be, [4]uint8{writer.MainVersion, MinorVersion, uint8(writer.K), 0})
	if writer.err != nil {
		return writer.err
	}

	writer.err = binary.Write(writer.w, be, writer.Flag)
	if writer.err != nil {
		return writer.err
	}
	return nil
}

// WriteKmer writes one Kmer.
func (writer *Writer) WriteKmer(mer []byte) error {
	writer.kcode, writer.err = NewKmerCode(mer)
	if writer.err != nil {
		return writer.err
	}
	return writer.Write(writer.kcode)
}

// Write writes one KmerCode.
func (writer *Writer) Write(kcode KmerCode) error {
	if writer.Header.K != kcode.K {
		writer.err = ErrKMismatch
		return writer.err
	}

	// lazily write header
	if !writer.wroteHeader {
		writer.err = writer.writeHeader()
		if writer.err != nil {
			return writer.err
		}
		writer.wroteHeader = true
	}

	if writer.compact {
		be.PutUint64(writer.buf, kcode.Code)
		writer.err = binary.Write(writer.w, be, writer.buf[8-writer.bufsize:])
	} else {
		writer.err = binary.Write(writer.w, be, kcode.Code)
	}
	if writer.err != nil {
		return writer.err
	}
	writer.size++
	return nil
}

// Flush is not used actually.
func (writer *Writer) Flush() error {
	// writer.err = binary.Write(writer.w, be, writer.size)
	// if writer.err != nil {
	// 	return writer.err
	// }
	return nil
}
