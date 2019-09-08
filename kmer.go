// Copyright © 2018-2019 Wei Shen <shenwei356@gmail.com>
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
//b
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

package unikmer

import (
	"bytes"
	"errors"
)

// ErrIllegalBase means that base beyond IUPAC symbols are  detected.
var ErrIllegalBase = errors.New("unikmer: illegal base")

// ErrKOverflow means K > 32.
var ErrKOverflow = errors.New("unikmer: k-mer size (1-32) overflow")

// ErrCodeOverflow means the encode interger is bigger than 4^k
var ErrCodeOverflow = errors.New("unikmer: code value overflow")

// slice is much faster than switch and map
var base2bit []uint64

// MaxCode is the maxinum interger for all Ks.
var MaxCode []uint64

func init() {
	MaxCode = make([]uint64, 33)
	for i := 1; i <= 32; i++ {
		MaxCode[i] = 1<<uint(i*2) - 1
	}

	base2bit = make([]uint64, 255)
	for i := range base2bit {
		base2bit[i] = 4
	}
	base2bit['A'] = 0
	base2bit['a'] = 0
	base2bit['N'] = 0
	base2bit['n'] = 0
	base2bit['M'] = 0
	base2bit['m'] = 0
	base2bit['V'] = 0
	base2bit['v'] = 0
	base2bit['H'] = 0
	base2bit['h'] = 0
	base2bit['R'] = 0
	base2bit['r'] = 0
	base2bit['D'] = 0
	base2bit['d'] = 0
	base2bit['W'] = 0
	base2bit['w'] = 0

	base2bit['C'] = 1
	base2bit['c'] = 1
	base2bit['S'] = 1
	base2bit['s'] = 1
	base2bit['B'] = 1
	base2bit['b'] = 1
	base2bit['Y'] = 1
	base2bit['y'] = 1

	base2bit['G'] = 2
	base2bit['g'] = 2
	base2bit['K'] = 2
	base2bit['k'] = 2

	base2bit['T'] = 3
	base2bit['t'] = 3
	base2bit['U'] = 3
	base2bit['u'] = 3
}

// Encode converts byte slice to bits.
//
// Codes:
//
// 	  A    0b00
// 	  C    0b01
// 	  G    0b10
// 	  T    0b11
//
// For degenerate bases, only the first base is kept.
//
//     M       AC     A
//     V       ACG    A
//     H       ACT    A
//     R       AG     A
//     D       AGT    A
//     W       AT     A
//     S       CG     C
//     B       CGT    C
//     Y       CT     C
//     K       GT     G
//     N       ACGT   A
//
func Encode(kmer []byte) (code uint64, err error) {
	if len(kmer) == 0 || len(kmer) > 32 {
		return 0, ErrKOverflow
	}

	var v uint64
	for _, b := range kmer {
		code <<= 2
		v = base2bit[b]
		if v > 3 {
			return code, ErrIllegalBase
		}
		code |= v
	}
	return code, nil
}

// ErrNotConsecutiveKmers means the two k-mers are not consecutive
var ErrNotConsecutiveKmers = errors.New("unikmer: not consecutive k-mers")

// MustEncodeFromFormerKmer encodes from former the k-mer,
// assuming the k-mer and leftKmer are both OK.
func MustEncodeFromFormerKmer(kmer []byte, leftKmer []byte, leftCode uint64) (uint64, error) {
	leftCode = leftCode & ((1 << (uint(len(kmer)-1) << 1)) - 1) << 2
	v := base2bit[kmer[len(kmer)-1]]
	if v > 3 {
		return leftCode, ErrIllegalBase
	}
	leftCode |= v
	return leftCode, nil
}

// EncodeFromFormerKmer encodes from the former k-mer, inspired by ntHash
func EncodeFromFormerKmer(kmer []byte, leftKmer []byte, leftCode uint64) (uint64, error) {
	if len(kmer) == 0 {
		return 0, ErrKOverflow
	}
	if len(kmer) != len(leftKmer) {
		return 0, ErrKMismatch
	}
	if !bytes.Equal(kmer[0:len(kmer)-1], leftKmer[1:len(leftKmer)]) {
		return 0, ErrNotConsecutiveKmers
	}
	return MustEncodeFromFormerKmer(kmer, leftKmer, leftCode)
}

// MustEncodeFromLatterKmer encodes from the latter k-mer,
// assuming the k-mer and rightKmer are both OK.
func MustEncodeFromLatterKmer(kmer []byte, rightKmer []byte, rightCode uint64) (uint64, error) {
	rightCode >>= 2
	v := base2bit[kmer[0]]
	if v > 3 {
		return rightCode, ErrIllegalBase
	}
	rightCode |= v << (uint(len(kmer)-1) << 1)
	return rightCode, nil
}

// EncodeFromLatterKmer encodes from the former k-mer.
func EncodeFromLatterKmer(kmer []byte, rightKmer []byte, rightCode uint64) (uint64, error) {
	if len(kmer) == 0 {
		return 0, ErrKOverflow
	}
	if len(kmer) != len(rightKmer) {
		return 0, ErrKMismatch
	}
	if !bytes.Equal(rightKmer[0:len(kmer)-1], kmer[1:len(rightKmer)]) {
		return 0, ErrNotConsecutiveKmers
	}
	return MustEncodeFromLatterKmer(kmer, rightKmer, rightCode)
}

// Reverse returns code of the reversed sequence.
func Reverse(code uint64, k int) (c uint64) {
	if k <= 0 || k > 32 {
		panic(ErrKOverflow)
	}
	for i := 0; i < k; i++ {
		c <<= 2
		c |= code & 3
		code >>= 2
	}
	return
}

// Complement returns code of complement sequence.
func Complement(code uint64, k int) (c uint64) {
	if k <= 0 || k > 32 {
		panic(ErrKOverflow)
	}
	for i := 0; i < k; i++ {
		c |= (code&3 ^ 3) << uint(i<<1)
		code >>= 2
	}
	return
}

// RevComp returns code of reverse complement sequence.
func RevComp(code uint64, k int) (c uint64) {
	if k <= 0 || k > 32 {
		panic(ErrKOverflow)
	}
	for i := 0; i < k; i++ {
		c <<= 2
		c |= code&3 ^ 3
		code >>= 2
	}
	return
}

// bit2base is for mapping bit to base.
var bit2base = [4]byte{'A', 'C', 'G', 'T'}

// bit2str is for output bits string
var bit2str = [4]string{"00", "01", "10", "11"}

// Decode converts the code to original seq
func Decode(code uint64, k int) []byte {
	if k <= 0 || k > 32 {
		panic(ErrKOverflow)
	}
	if code > MaxCode[k] {
		panic(ErrCodeOverflow)
	}
	kmer := make([]byte, k)
	for i := 0; i < k; i++ {
		kmer[k-1-i] = bit2base[code&3]
		code >>= 2
	}
	return kmer
}

// KmerCode is a struct representing a k-mer in 64-bits.
type KmerCode struct {
	Code uint64
	K    int
}

// NewKmerCode returns a new KmerCode struct from byte slice.
func NewKmerCode(kmer []byte) (KmerCode, error) {
	code, err := Encode(kmer)
	if err != nil {
		return KmerCode{}, err
	}
	return KmerCode{code, len(kmer)}, err
}

// NewKmerCodeFromFormerOne computes KmerCode from the Former consecutive k-mer.
func NewKmerCodeFromFormerOne(kmer []byte, leftKmer []byte, preKcode KmerCode) (KmerCode, error) {
	code, err := EncodeFromFormerKmer(kmer, leftKmer, preKcode.Code)
	if err != nil {
		return KmerCode{}, err
	}
	return KmerCode{code, len(kmer)}, err
}

// NewKmerCodeMustFromFormerOne computes KmerCode from the Former consecutive k-mer,
// assuming the k-mer and leftKmer are both OK.
func NewKmerCodeMustFromFormerOne(kmer []byte, leftKmer []byte, preKcode KmerCode) (KmerCode, error) {
	code, err := MustEncodeFromFormerKmer(kmer, leftKmer, preKcode.Code)
	if err != nil {
		return KmerCode{}, err
	}
	return KmerCode{code, len(kmer)}, err
}

// Equal checks wether two KmerCodes are the same.
func (kcode KmerCode) Equal(kcode2 KmerCode) bool {
	return kcode.K == kcode2.K && kcode.Code == kcode2.Code
}

// Rev returns KmerCode of the reverse sequence.
func (kcode KmerCode) Rev() KmerCode {
	return KmerCode{Reverse(kcode.Code, kcode.K), kcode.K}
}

// Comp returns KmerCode of the complement sequence.
func (kcode KmerCode) Comp() KmerCode {
	return KmerCode{Complement(kcode.Code, kcode.K), kcode.K}
}

// RevComp returns KmerCode of the reverse complement sequence.
func (kcode KmerCode) RevComp() KmerCode {
	return KmerCode{RevComp(kcode.Code, kcode.K), kcode.K}
}

// Canonical returns its canonical kmer
func (kcode KmerCode) Canonical() KmerCode {
	rcKcode := kcode.RevComp()
	if rcKcode.Code < kcode.Code {
		return rcKcode
	}
	return kcode
}

// Bytes returns k-mer in []byte.
func (kcode KmerCode) Bytes() []byte {
	return Decode(kcode.Code, kcode.K)
}

// String returns k-mer in string
func (kcode KmerCode) String() string {
	return string(Decode(kcode.Code, kcode.K))
}

// BitsString returns code to string
func (kcode KmerCode) BitsString() string {
	var buf bytes.Buffer
	for _, b := range Decode(kcode.Code, kcode.K) {
		buf.WriteString(bit2str[base2bit[b]])
	}
	return buf.String()
}
