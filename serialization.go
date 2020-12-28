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

package unikmer

import (
	"encoding/binary"
	"errors"
	"fmt"
	"io"
)

// MainVersion is the main version number.
const MainVersion uint8 = 5

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

// ErrDescTooLong means length of description two long
var ErrDescTooLong = errors.New("unikmer: description too long, 128 bytes at most")

// ErrCallOrder means WriteTaxid/ReadTaxid should be called after WriteCode/ReadCode
var ErrCallOrder = errors.New("unikmer: WriteTaxid/ReadTaxid should be called after WriteCode/ReadCode")

// ErrCallLate means SetMaxTaxid/SetGlobalTaxid should be called before writing KmerCode/code/taxid
var ErrCallLate = errors.New("unikmer: SetMaxTaxid/SetGlobalTaxid should be called before writing KmerCode/code/taxid")

// ErrCallReadWriteTaxid means flag UnikIncludeTaxID is off, but you call ReadTaxid/WriteTaxid
var ErrCallReadWriteTaxid = errors.New("unikmer: can not call ReadTaxid/WriteTaxid when flag UnikIncludeTaxID is off")

// ErrInvalidTaxid means zero given for a taxid.
var ErrInvalidTaxid = errors.New("unikmer: invalid taxid, 0 not allowed")

// ErrVersionMismatch means version mismatch between files and program
var ErrVersionMismatch = errors.New("unikmer: version mismatch")

var be = binary.BigEndian

var descMaxLen = 1024
var conservedDataLen = 64

// Header contains metadata
type Header struct {
	MainVersion  uint8
	MinorVersion uint8
	K            int
	Flag         uint32
	Number       int64  // -1 for unknown
	globalTaxid  uint32 // universal taxid, 0 for no record
	maxTaxid     uint32
	Description  []byte // let's limit it to 128 Bytes
	Scale        uint32 // scale of down-sampling
	MaxHash      uint64 // max hash for scaling/down-sampling
}

const (
	// UnikCompact means k-mers are serialized in fix-length (n = int((K + 3) / 4) ) of byte array.
	UnikCompact = 1 << iota
	// UnikCanonical means only canonical k-mers kept.
	UnikCanonical
	// UnikSorted means k-mers are sorted
	UnikSorted // when sorted, the serialization structure is very different
	// UnikIncludeTaxID means a k-mer are followed it's LCA taxid
	UnikIncludeTaxID

	// UnikHashed means ntHash value are saved as code.
	UnikHashed
	// UnikScaled means only hashes smaller than or equal to max_hash are saved.
	UnikScaled
)

func (h Header) String() string {
	return fmt.Sprintf("unikmer binary k-mer data file v%d.%d with K=%d and Flag=%d",
		h.MainVersion, h.MinorVersion, h.K, h.Flag)
}

// Reader is for reading KmerCode.
type Reader struct {
	Header
	r io.Reader

	buf []byte

	compact bool // saving KmerCode in variable-length byte array.
	bufsize int

	sorted  bool
	hasPrev bool
	prev    uint64
	buf2    []byte
	offset  uint64

	includeTaxid  bool
	bufTaxid      []byte
	taxidByteLen  int
	prevTaxid     uint32 // buffered taxid
	hasPrevTaxid  bool
	justReadACode bool
	lastRecord    bool

	hashValue bool
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

// IsSorted tells if the k-mers in file sorted
func (reader *Reader) IsSorted() bool {
	return reader.Flag&UnikSorted > 0
}

// IsCanonical tells if the only canonical k-mers stored
func (reader *Reader) IsCanonical() bool {
	return reader.Flag&UnikCanonical > 0
}

// IsCompact tells if the k-mers are stored in a compact format
func (reader *Reader) IsCompact() bool {
	return reader.Flag&UnikCompact > 0
}

// IsIncludeTaxid tells if every k-mer is followed by its taxid
func (reader *Reader) IsIncludeTaxid() bool {
	return reader.Flag&UnikIncludeTaxID > 0
}

// IsHashed tells if ntHash values are saved.
func (reader *Reader) IsHashed() bool {
	return reader.Flag&UnikHashed > 0
}

// IsScaled tells if hashes is scaled
func (reader *Reader) IsScaled() bool {
	return reader.Flag&UnikHashed > 0 && reader.Flag&UnikScaled > 0
}

// HasGlobalTaxid means the file has a global taxid
func (reader *Reader) HasGlobalTaxid() bool {
	return reader.globalTaxid > 0
}

// HasTaxidInfo means the binary file contains global taxid or taxids for all k-mers
func (reader *Reader) HasTaxidInfo() bool {
	return reader.IsIncludeTaxid() || reader.HasGlobalTaxid()
}

// GetGlobalTaxid returns the global taxid
func (reader *Reader) GetGlobalTaxid() uint32 {
	return reader.globalTaxid
}

// GetTaxidBytesLength returns number of byte to store a taxid
func (reader *Reader) GetTaxidBytesLength() int {
	return reader.taxidByteLen
}

// GetScale returns the scale of down-sampling
func (reader *Reader) GetScale() uint32 {
	if reader.Scale == 0 {
		return uint32(1)
	}
	return reader.Scale
}

// GetMaxHash returns the max hash for scaling.
func (reader *Reader) GetMaxHash() uint64 {
	if reader.MaxHash == 0 {
		return ^uint64(0)
	}
	return reader.MaxHash
}

func (reader *Reader) readHeader() (err error) {
	buf := make([]byte, 56)
	r := reader.r

	// check Magic number
	_, err = io.ReadFull(r, buf[:8])
	if err != nil {
		return err
	}
	same := true
	for i := 0; i < 8; i++ {
		if Magic[i] != buf[i] {
			same = false
			break
		}
	}
	if !same {
		return ErrInvalidFileFormat
	}

	// read metadata
	_, err = io.ReadFull(r, buf[:4])
	if err != nil {
		return err
	}
	// check compatibility？
	if (buf[0] == 0 && buf[1] == 0) ||
		MainVersion != buf[0] {
		return ErrVersionMismatch
	}
	reader.MainVersion = buf[0]
	reader.MinorVersion = buf[1]

	reader.K = int(buf[2])

	_, err = io.ReadFull(r, buf[:4])
	if err != nil {
		return err
	}
	reader.Flag = be.Uint32(buf[:4])

	reader.buf = make([]byte, 8)

	if reader.IsCompact() {
		reader.compact = true
		reader.bufsize = int((reader.K + 3) / 4)
	}
	if reader.IsSorted() {
		reader.sorted = true
		reader.buf2 = make([]byte, 17)
	}
	if reader.IsIncludeTaxid() {
		reader.includeTaxid = true
		reader.bufTaxid = make([]byte, 4)
	}

	// number
	_, err = io.ReadFull(r, buf[:8])
	if err != nil {
		return err
	}
	reader.Number = int64(be.Uint64(buf[:8]))

	// taxid
	_, err = io.ReadFull(r, buf[:4])
	if err != nil {
		return err
	}
	reader.globalTaxid = be.Uint32(buf[:4])

	// taxid byte length
	_, err = io.ReadFull(r, buf[:1])
	if err != nil {
		return err
	}
	buf[1] = 0
	reader.taxidByteLen = int(be.Uint16(buf[:2]))

	// length of description
	var lenDesc uint16
	_, err = io.ReadFull(r, buf[:2])
	if err != nil {
		return err
	}
	lenDesc = be.Uint16(buf[:2])

	desc := make([]byte, lenDesc)
	_, err = io.ReadFull(r, desc)
	if err != nil {
		return err
	}
	reader.Description = desc

	// scale
	_, err = io.ReadFull(r, buf[:4])
	if err != nil {
		return err
	}
	reader.Scale = be.Uint32(buf[:4])

	// max hash
	_, err = io.ReadFull(r, buf[:8])
	if err != nil {
		return err
	}
	reader.MaxHash = be.Uint64(buf[:8])

	reserved := make([]byte, conservedDataLen)
	_, err = io.ReadFull(r, reserved)
	if err != nil {
		return err
	}

	return nil
}

// Read reads one KmerCode.
func (reader *Reader) Read() (KmerCode, error) {
	code, err := reader.ReadCode()
	return KmerCode{Code: code, K: reader.K}, err
}

// ReadWithTaxid reads a KmerCode, also return taxid if having.
func (reader *Reader) ReadWithTaxid() (KmerCode, uint32, error) {
	code, taxid, err := reader.ReadCodeWithTaxid()
	return KmerCode{Code: code, K: reader.K}, taxid, err
}

// ReadCodeWithTaxid reads a code, also return taxid if having.
func (reader *Reader) ReadCodeWithTaxid() (code uint64, taxid uint32, err error) {
	code, err = reader.ReadCode()
	if err != nil {
		return 0, 0, err
	}
	if reader.includeTaxid {
		taxid, err = reader.ReadTaxid()
		if err != nil {
			return 0, 0, err
		}
	} else {
		taxid = reader.globalTaxid
	}
	return code, taxid, err
}

// ReadTaxid reads on taxid
func (reader *Reader) ReadTaxid() (taxid uint32, err error) {
	if !reader.includeTaxid {
		return 0, ErrCallReadWriteTaxid
	}

	if !reader.justReadACode {
		return 0, ErrCallOrder
	}

	if reader.sorted {
		if reader.lastRecord {
			_, err = io.ReadFull(reader.r, reader.bufTaxid)
			if err != nil {
				return 0, err
			}
			reader.hasPrevTaxid = false
			reader.justReadACode = false
			return be.Uint32(reader.bufTaxid), nil
		}

		if reader.hasPrevTaxid {
			c := reader.prevTaxid
			reader.hasPrevTaxid = false
			reader.justReadACode = false
			return c, nil
		}

		_, err = io.ReadFull(reader.r, reader.bufTaxid[4-reader.taxidByteLen:])
		if err != nil {
			return 0, err
		}
		taxid = be.Uint32(reader.bufTaxid)

		_, err = io.ReadFull(reader.r, reader.bufTaxid[4-reader.taxidByteLen:])
		if err != nil {
			return 0, err
		}

		reader.prevTaxid = be.Uint32(reader.bufTaxid)
		reader.hasPrevTaxid = true
		return taxid, nil
	} else if reader.compact {
		_, err = io.ReadFull(reader.r, reader.bufTaxid[4-reader.taxidByteLen:])
	} else {
		_, err = io.ReadFull(reader.r, reader.bufTaxid)
	}
	if err != nil {
		return 0, err
	}

	reader.justReadACode = false
	return be.Uint32(reader.bufTaxid), nil
}

// ReadCode reads one code.
func (reader *Reader) ReadCode() (uint64, error) {
	var err error
	if reader.sorted {
		if reader.hasPrev {
			c := reader.prev
			// reader.prev = 0
			reader.hasPrev = false
			reader.justReadACode = true
			return c, nil
		}

		buf2 := reader.buf2
		r := reader.r

		// read control byte
		var nReaded int
		nReaded, err = io.ReadFull(r, buf2[0:1])
		if err != nil {
			return 0, err
		}

		ctrlByte := buf2[0]
		if ctrlByte&128 > 0 { // last one
			nReaded, err = io.ReadFull(r, buf2[0:8])
			if err != nil {
				return 0, err
			}
			reader.lastRecord = true
			reader.justReadACode = true
			return be.Uint64(buf2[0:8]), nil
		}

		// parse control byte
		encodedBytes := ctrlByte2ByteLengths[ctrlByte]
		nEncodedBytes := int(encodedBytes[0] + encodedBytes[1])

		// read encoded bytes
		nReaded, err = io.ReadFull(r, buf2[0:nEncodedBytes])
		if err != nil {
			return 0, err
		}
		if nReaded < nEncodedBytes {
			return 0, ErrBrokenFile
		}

		decodedVals, nDecoded := Uint64s(ctrlByte, buf2[0:nEncodedBytes])
		if nDecoded == 0 {
			return 0, ErrBrokenFile
		}

		code := decodedVals[0] + reader.offset
		reader.prev = code + decodedVals[1]
		reader.hasPrev = true

		reader.offset = code + decodedVals[1]

		reader.justReadACode = true
		return code, nil
	} else if reader.compact {
		_, err = io.ReadFull(reader.r, reader.buf[8-reader.bufsize:])
	} else {
		_, err = io.ReadFull(reader.r, reader.buf)
	}
	if err != nil {
		return 0, err
	}

	reader.justReadACode = true
	return be.Uint64(reader.buf), nil
}

// Writer writes KmerCode.
type Writer struct {
	Header
	w           io.Writer
	wroteHeader bool

	buf []byte

	// saving KmerCode in compact fixlength byte array.
	compact bool
	bufsize int

	// sortred mode
	sorted       bool
	offset       uint64
	prev         uint64 // buffered code
	hasPrev      bool
	buf2         []byte
	buf3         []byte
	ctrlByte     byte
	nEncodedByte int

	// for taxid
	includeTaxid     bool
	bufTaxid         []byte
	justWrittenACode bool
	taxidByteLen     int
	prevTaxid        uint32 // buffered taxid
	hasPrevTaxid     bool
}

// NewWriter creates a Writer.
func NewWriter(w io.Writer, k int, flag uint32) (*Writer, error) {
	if k == 0 {
		return nil, ErrKOverflow
	}

	writer := &Writer{
		Header: Header{MainVersion: MainVersion, MinorVersion: MinorVersion, K: k, Flag: flag, Number: -1},
		w:      w,
	}

	// prevent wrong use of compact
	if writer.Flag&UnikCompact > 0 && (writer.Flag&UnikSorted > 0 || writer.Flag&UnikHashed > 0) {
		writer.Flag ^= UnikCompact
	}

	writer.buf = make([]byte, 8)
	if writer.Flag&UnikCompact > 0 &&
		writer.Flag&UnikSorted == 0 &&
		writer.Flag&UnikHashed == 0 {

		writer.compact = true
		writer.bufsize = int(k+3) / 4
	} else if writer.Flag&UnikSorted > 0 {
		writer.sorted = true
		writer.buf2 = make([]byte, 16)
		writer.buf3 = make([]byte, 32)
	}
	if writer.Flag&UnikIncludeTaxID > 0 {
		writer.includeTaxid = true
		writer.bufTaxid = make([]byte, 4)
	}

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
	err = binary.Write(w, be, [4]uint8{writer.MainVersion, MinorVersion, uint8(writer.K), 0})
	if err != nil {
		return err
	}

	// 4 bytes flags
	err = binary.Write(w, be, writer.Flag)
	if err != nil {
		return err
	}

	// 8 bytes number
	err = binary.Write(w, be, writer.Number)
	if err != nil {
		return err
	}

	// 4 bytes taxid
	err = binary.Write(w, be, writer.globalTaxid)
	if err != nil {
		return err
	}

	// 1 byte taxid bytes len
	if writer.maxTaxid <= 0 {
		writer.taxidByteLen = 4
	} else {
		writer.taxidByteLen = int(byteLength(uint64(writer.maxTaxid)))
	}
	err = binary.Write(w, be, uint8(writer.taxidByteLen))
	if err != nil {
		return err
	}

	// description length (2 byte)s and data (128 bytes)
	lenDesc := len(writer.Description)
	if lenDesc > descMaxLen {
		return ErrDescTooLong
	}
	err = binary.Write(w, be, uint16(lenDesc))
	if err != nil {
		return err
	}
	err = binary.Write(w, be, writer.Description)
	if err != nil {
		return err
	}

	// scale
	err = binary.Write(w, be, writer.Scale)
	if err != nil {
		return err
	}

	// max hash
	err = binary.Write(w, be, writer.MaxHash)
	if err != nil {
		return err
	}

	// reserved 24 bytes
	reserved := make([]byte, conservedDataLen)
	err = binary.Write(w, be, reserved)
	if err != nil {
		return err
	}

	// header has 192 bytes

	writer.wroteHeader = true
	return nil
}

// SetGlobalTaxid sets the global taxid
func (writer *Writer) SetGlobalTaxid(taxid uint32) error {
	if writer.wroteHeader {
		return ErrCallLate
	}
	writer.globalTaxid = taxid
	return nil
}

// SetMaxTaxid set the maxtaxid
func (writer *Writer) SetMaxTaxid(taxid uint32) error {
	if writer.wroteHeader {
		return ErrCallLate
	}
	writer.maxTaxid = taxid
	return nil
}

// SetScale set the scale
func (writer *Writer) SetScale(scale uint32) error {
	if writer.wroteHeader {
		return ErrCallLate
	}
	if writer.Flag&UnikHashed == 0 {
		writer.Flag |= UnikHashed
	}
	if writer.Flag&UnikScaled == 0 {
		writer.Flag |= UnikScaled
	}
	writer.Scale = scale
	return nil
}

// SetMaxHash set the max hash
func (writer *Writer) SetMaxHash(maxHash uint64) error {
	if writer.wroteHeader {
		return ErrCallLate
	}
	if writer.Flag&UnikHashed == 0 {
		writer.Flag += UnikHashed
	}
	if writer.Flag&UnikScaled == 0 {
		writer.Flag += UnikScaled
	}
	writer.MaxHash = maxHash
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

// WriteKmerWithTaxid writes one k-mer and its taxid
func (writer *Writer) WriteKmerWithTaxid(mer []byte, taxid uint32) error {
	err := writer.WriteKmer(mer)
	if err != nil {
		return nil
	}
	return writer.WriteTaxid(taxid)
}

// Write writes one KmerCode.
func (writer *Writer) Write(kcode KmerCode) (err error) {
	if writer.K != kcode.K {
		return ErrKMismatch
	}
	return writer.WriteCode(kcode.Code)
}

// WriteWithTaxid writes one KmerCode and its taxid.
// If UnikIncludeTaxID is off, taxid will not be written.
func (writer *Writer) WriteWithTaxid(kcode KmerCode, taxid uint32) (err error) {
	err = writer.Write(kcode)
	if err != nil {
		return nil
	}
	return writer.WriteTaxid(taxid)
}

// WriteCodeWithTaxid writes a code and its taxid.
// If UnikIncludeTaxID is off, taxid will not be written.
func (writer *Writer) WriteCodeWithTaxid(code uint64, taxid uint32) (err error) {
	err = writer.WriteCode(code)
	if err != nil {
		return nil
	}
	if !writer.includeTaxid { // if no taxid, just return.
		return nil
	}
	return writer.WriteTaxid(taxid)
}

// WriteTaxid appends taxid to the code
func (writer *Writer) WriteTaxid(taxid uint32) (err error) {
	if !writer.includeTaxid {
		return ErrCallReadWriteTaxid
	}

	if !writer.justWrittenACode {
		return ErrCallOrder
	}

	if writer.sorted {
		if !writer.hasPrevTaxid { // write it later
			writer.prevTaxid = taxid
			writer.hasPrevTaxid = true
			writer.justWrittenACode = false
			return nil
		}
		be.PutUint32(writer.bufTaxid, writer.prevTaxid)
		_, err = writer.w.Write(writer.bufTaxid[4-writer.taxidByteLen:])
		// fmt.Printf("write taxid: %d, %d\n", writer.prevTaxid, writer.bufTaxid[4-writer.taxidByteLen:])

		be.PutUint32(writer.bufTaxid, taxid)
		_, err = writer.w.Write(writer.bufTaxid[4-writer.taxidByteLen:])
		writer.hasPrevTaxid = false
	} else if writer.compact {
		be.PutUint32(writer.bufTaxid, taxid)
		_, err = writer.w.Write(writer.bufTaxid[4-writer.taxidByteLen:])
	} else {
		be.PutUint32(writer.bufTaxid, taxid)
		_, err = writer.w.Write(writer.bufTaxid)
	}

	writer.justWrittenACode = false
	return nil
}

// WriteCode writes one code
func (writer *Writer) WriteCode(code uint64) (err error) {
	// lazily write header
	if !writer.wroteHeader {
		err = writer.WriteHeader()
		if err != nil {
			return err
		}
		writer.wroteHeader = true
	}

	if writer.sorted {
		if !writer.hasPrev { // write it later
			writer.prev = code
			writer.hasPrev = true
			writer.justWrittenACode = true
			return nil
		}

		writer.ctrlByte, writer.nEncodedByte = PutUint64s(writer.buf2, writer.prev-writer.offset, code-writer.prev)

		writer.buf3[0] = writer.ctrlByte
		copy(writer.buf3[1:writer.nEncodedByte+1], writer.buf2[0:writer.nEncodedByte])
		_, err = writer.w.Write(writer.buf3[0 : writer.nEncodedByte+1])

		writer.offset = code
		// writer.prev = 0
		writer.hasPrev = false
	} else if writer.compact {
		be.PutUint64(writer.buf, code)
		_, err = writer.w.Write(writer.buf[8-writer.bufsize:])
	} else {
		be.PutUint64(writer.buf, code)
		_, err = writer.w.Write(writer.buf)
	}

	if err != nil {
		return err
	}
	writer.justWrittenACode = true
	return nil
}

// Flush write the last k-mer
func (writer *Writer) Flush() (err error) {
	if !writer.wroteHeader {
		writer.Number = 0
		writer.WriteHeader()
	}
	if !writer.sorted || !writer.hasPrev {
		return nil
	}

	// write last k-mer
	err = binary.Write(writer.w, be, uint8(128))
	if err != nil {
		return err
	}
	err = binary.Write(writer.w, be, writer.prev) // last code
	if err != nil {
		return err
	}
	if writer.includeTaxid && writer.hasPrevTaxid { // last taxid
		err = binary.Write(writer.w, be, writer.prevTaxid)
		if err != nil {
			return err
		}
	}

	writer.hasPrev = false
	writer.hasPrevTaxid = false
	return nil
}
