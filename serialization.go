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
var ErrBrokenFile = errors.New("unikmer: broken file")

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
	// UNIK_COMPACT means Kmers are serialized in fix-length (n = int((K + 3) / 4) ) of byte array.
	UNIK_COMPACT = 1 << iota
	// UNIK_CANONICAL means only canonical Kmers kept.
	UNIK_CANONICAL
	// UNIK_SORTED means Kmers are sorted
	UNIK_SORTED // when sorted, the serialization structure is very different
)

func (h Header) String() string {
	return fmt.Sprintf("unikmer binary k-mer data file v%d.%d with K=%d and Flag=%d",
		h.MainVersion, h.MinorVersion, h.K, h.Flag)
}

// Reader is for reading KmerCode.
type Reader struct {
	Header
	r   io.Reader
	buf []byte

	compact bool // saving KmerCode in variable-length byte array.
	bufsize int

	sorted bool
	prev   *KmerCode
	buf2   []byte
	offset uint64
}

// NewReader returns a Reader.
func NewReader(r io.Reader) (reader *Reader, err error) {
	reader = &Reader{r: r}
	err = reader.readHeader()
	if err != nil {
		return nil, err
	}
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
		return ErrInvalidFileFormat
	}

	// read metadata
	var meta [4]uint8
	err = binary.Read(r, be, &meta)
	if err != nil {
		return err
	}
	// check compatibility？
	if (meta[0] == 0 && meta[1] == 0) ||
		MainVersion != meta[0] {
		return fmt.Errorf("unikmer: .unik format compatibility error, please recreate with newest version")
	}
	reader.MainVersion = meta[0]
	reader.MinorVersion = meta[1]

	reader.K = int(meta[2])

	err = binary.Read(r, be, &reader.Flag)
	if err != nil {
		return err
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

	err = binary.Read(r, be, &reader.Number)
	if err != nil {
		return err
	}
	return nil
}

// Read reads one KmerCode.
func (reader *Reader) Read() (KmerCode, error) {
	var err error
	if reader.sorted {
		if reader.prev != nil {
			c := *reader.prev
			reader.prev = nil
			return c, nil
		}

		buf2 := reader.buf2
		r := reader.r

		// read control byte
		var nReaded int
		nReaded, err = io.ReadFull(r, buf2[0:1])
		if err != nil {
			return KmerCode{}, err
		}

		ctrlByte := buf2[0]
		if ctrlByte&128 > 0 { // last one
			nReaded, err = io.ReadFull(r, buf2[0:8])
			if err != nil {
				return KmerCode{}, err
			}
			return KmerCode{Code: be.Uint64(buf2[0:8]), K: reader.K}, nil
		}

		// parse control byte
		encodedBytes := ctrlByte2ByteLengths[ctrlByte]
		nEncodedBytes := int(encodedBytes[0] + encodedBytes[1])

		// read encoded bytes
		nReaded, err = io.ReadFull(r, buf2[0:nEncodedBytes])
		if err != nil {
			return KmerCode{}, err
		}
		if nReaded < nEncodedBytes {
			return KmerCode{}, ErrBrokenFile
		}

		decodedVals, nDecoded := Uint64s(ctrlByte, buf2[0:nEncodedBytes])
		if nDecoded == 0 {
			return KmerCode{}, ErrBrokenFile
		}

		code := decodedVals[0] + reader.offset
		reader.prev = &KmerCode{Code: code + decodedVals[1], K: reader.K}
		reader.offset = code + decodedVals[1]

		return KmerCode{Code: code, K: reader.K}, nil
	} else if reader.compact {
		_, err = io.ReadFull(reader.r, reader.buf[8-reader.bufsize:])
	} else {
		_, err = io.ReadFull(reader.r, reader.buf)
	}
	if err != nil {
		return KmerCode{}, err
	}

	return KmerCode{Code: be.Uint64(reader.buf), K: reader.K}, nil
}

// Writer writes KmerCode.
type Writer struct {
	Header
	w           io.Writer
	wroteHeader bool

	compact bool // saving KmerCode in compact fixlength byte array.
	buf     []byte
	bufsize int

	sorted bool //
	prev   *KmerCode
	buf2   []byte
	offset uint64
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
func (writer *Writer) WriteHeader() (err error) {
	if writer.wroteHeader {
		return nil
	}
	w := writer.w
	// write magic number
	err = binary.Write(w, be, Magic)
	if err != nil {
		return err
	}

	err = binary.Write(w, be, [4]uint8{writer.MainVersion, MinorVersion, uint8(writer.K), 0})
	if err != nil {
		return err
	}

	err = binary.Write(w, be, writer.Flag)
	if err != nil {
		return err
	}

	err = binary.Write(w, be, writer.Number)
	if err != nil {
		return err
	}

	writer.wroteHeader = true
	return nil
}

// WriteKmer writes one k-mer.
func (writer *Writer) WriteKmer(mer []byte) error {
	kcode, err := NewKmerCode(mer)
	if err != nil {
		return err
	}
	return writer.Write(kcode)
}

// Write writes one KmerCode.
func (writer *Writer) Write(kcode KmerCode) (err error) {
	if writer.K != kcode.K {
		return ErrKMismatch
	}

	// lazily write header
	if !writer.wroteHeader {
		err = writer.WriteHeader()
		if err != nil {
			return err
		}
		writer.wroteHeader = true
	}

	if writer.sorted {
		if writer.prev == nil { // write it later
			writer.prev = &kcode
			return nil
		}

		ctrlByte, nEncodedByte := PutUint64s(writer.buf2, writer.prev.Code-writer.offset, kcode.Code-writer.prev.Code)

		writer.offset = kcode.Code
		writer.prev = nil

		err = binary.Write(writer.w, be, ctrlByte)
		err = binary.Write(writer.w, be, writer.buf2[0:nEncodedByte])
	} else if writer.compact {
		be.PutUint64(writer.buf, kcode.Code)
		err = binary.Write(writer.w, be, writer.buf[8-writer.bufsize:])
	} else {
		err = binary.Write(writer.w, be, kcode.Code)
	}

	if err != nil {
		return err
	}
	return nil
}

// Flush write the last k-mer
func (writer *Writer) Flush() (err error) {
	if !writer.sorted || writer.prev == nil {
		return nil
	}
	// write last k-mer
	err = binary.Write(writer.w, be, uint8(128))
	err = binary.Write(writer.w, be, writer.prev.Code)
	if err != nil {
		return err
	}
	writer.prev = nil
	return nil
}
