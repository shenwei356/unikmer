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
const MainVersion uint8 = 2

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
	Number       int64 // -1 for unknown
}

const (
	UNIK_COMPACT = 1 << iota
	UNIK_CANONICAL
	UNIK_SORTED // when sorted, the serialization structure is very different
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

	buf []byte

	compact bool // saving KmerCode in variable-length byte array.
	bufsize int

	sorted        bool
	prev          *KmerCode
	buf2          []byte
	offset        uint64
	ctrlByte      byte
	nReaded       int
	nEncodedBytes int
	encodedBytes  [2]uint8
	decodedVals   [2]uint64
	nDecoded      int
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
		return fmt.Errorf("unikmer: .unik format compatibility error, please recreate with newest version")
	}
	reader.MainVersion = meta[0]
	reader.MinorVersion = meta[1]

	reader.K = int(meta[2])

	reader.err = binary.Read(reader.r, be, &reader.Flag)
	if reader.err != nil {
		return reader.err
	}

	reader.buf = make([]byte, 8)

	if reader.Flag&UNIK_COMPACT > 0 {
		reader.compact = true
		reader.bufsize = int((reader.K + 3) / 4)
	}
	if reader.Flag&UNIK_SORTED > 0 {
		reader.sorted = true
		reader.buf2 = make([]byte, 17)
	}

	reader.err = binary.Read(reader.r, be, &reader.Number)
	if reader.err != nil {
		return reader.err
	}
	return nil
}

// Read reads one KmerCode.
func (reader *Reader) Read() (KmerCode, error) {
	if reader.sorted {
		if reader.prev != nil {
			c := *reader.prev
			reader.prev = nil
			return c, nil
		}
		// read control byte
		reader.nReaded, reader.err = io.ReadFull(reader.r, reader.buf2[0:1])
		if reader.err != nil {
			return KmerCode{}, reader.err
		}

		reader.ctrlByte = reader.buf2[0]
		if reader.ctrlByte&128 > 0 { // last one
			reader.nReaded, reader.err = io.ReadFull(reader.r, reader.buf2[0:8])
			if reader.err != nil {
				return KmerCode{}, reader.err
			}
			reader.code = be.Uint64(reader.buf2[0:8])
			return KmerCode{Code: reader.code, K: reader.Header.K}, nil
		}

		// parse control byte
		reader.encodedBytes = ctrlByte2ByteLengths[reader.ctrlByte]
		reader.nEncodedBytes = int(reader.encodedBytes[0] + reader.encodedBytes[1])

		// read encoded bytes
		reader.nReaded, reader.err = io.ReadFull(reader.r, reader.buf2[0:reader.nEncodedBytes])
		if reader.err != nil {
			return KmerCode{}, reader.err
		}
		if reader.nReaded < reader.nEncodedBytes {
			return KmerCode{}, fmt.Errorf("unikmer: unexpected EOF")
		}

		reader.decodedVals, reader.nDecoded = Uint64s(reader.ctrlByte, reader.buf2[0:reader.nEncodedBytes])
		if reader.nDecoded == 0 {
			return KmerCode{Code: reader.code, K: reader.Header.K}, fmt.Errorf("unikmer: broken binary file")
		}

		reader.code = reader.decodedVals[0] + reader.offset
		reader.prev = &KmerCode{Code: reader.code + reader.decodedVals[1], K: reader.Header.K}
		reader.offset = reader.code + reader.decodedVals[1]

		return KmerCode{Code: reader.code, K: reader.Header.K}, nil
	} else if reader.compact {
		_, reader.err = io.ReadFull(reader.r, reader.buf[8-reader.bufsize:])
	} else {
		_, reader.err = io.ReadFull(reader.r, reader.buf)
	}
	if reader.err != nil {
		return KmerCode{}, reader.err
	}
	reader.code = be.Uint64(reader.buf)

	return KmerCode{Code: reader.code, K: reader.Header.K}, nil
}

// Writer writes KmerCode.
type Writer struct {
	Header
	w           io.Writer
	kcode       KmerCode
	wroteHeader bool
	err         error

	compact bool // saving KmerCode in compact fixlength byte array.
	buf     []byte
	bufsize int

	sorted       bool //
	prev         *KmerCode
	buf2         []byte
	offset       uint64
	ctrlByte     byte
	nEncodedByte int
}

// NewWriter creates a Writer.
func NewWriter(w io.Writer, k int, flag uint32) (*Writer, error) {
	if k == 0 || k > 32 {
		return nil, ErrKOverflow
	}

	writer := &Writer{
		Header: Header{MainVersion: MainVersion, MinorVersion: MinorVersion, K: k, Flag: flag, Number: -1},
		w:      w,
	}

	writer.buf = make([]byte, 8)
	if writer.Flag&UNIK_COMPACT > 0 {
		writer.compact = true
		writer.bufsize = int(k+3) / 4
	}
	if writer.Flag&UNIK_SORTED > 0 {
		writer.sorted = true
		writer.buf2 = make([]byte, 16)
	}
	return writer, nil
}

// WriteHeader writes file header
func (writer *Writer) WriteHeader() error {
	if writer.wroteHeader {
		return nil
	}
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

	writer.err = binary.Write(writer.w, be, writer.Number)
	if writer.err != nil {
		return writer.err
	}

	writer.wroteHeader = true
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
		writer.err = writer.WriteHeader()
		if writer.err != nil {
			return writer.err
		}
		writer.wroteHeader = true
	}

	if writer.sorted {
		if writer.prev == nil { // write it later
			writer.prev = &kcode
			return nil
		}

		writer.ctrlByte, writer.nEncodedByte = PutUint64s(writer.buf2, writer.prev.Code-writer.offset, kcode.Code-writer.prev.Code)

		writer.offset = kcode.Code
		writer.prev = nil

		writer.err = binary.Write(writer.w, be, writer.ctrlByte)
		writer.err = binary.Write(writer.w, be, writer.buf2[0:writer.nEncodedByte])
	} else if writer.compact {
		be.PutUint64(writer.buf, kcode.Code)
		writer.err = binary.Write(writer.w, be, writer.buf[8-writer.bufsize:])
	} else {
		writer.err = binary.Write(writer.w, be, kcode.Code)
	}

	if writer.err != nil {
		return writer.err
	}
	return nil
}

// Flush write the last Kmer
func (writer *Writer) Flush() error {
	if !writer.sorted || writer.prev == nil {
		return nil
	}
	// write last kmer
	writer.err = binary.Write(writer.w, be, uint8(128))
	writer.err = binary.Write(writer.w, be, writer.prev.Code)
	if writer.err != nil {
		return writer.err
	}
	writer.prev = nil
	return nil
}
