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
	"bytes"
	"errors"
)

// ErrIllegalBase means that base beyond IUPAC symbols are  detected.
var ErrIllegalBase = errors.New("unikmer: illegal base")

// ErrKOverflow means K > 32.
var ErrKOverflow = errors.New("unikmer: K (1-32) overflow")

// Encode converts byte slice to bits.
//
// Codes:
//
// 	  A    00
// 	  C    01
// 	  G    10
// 	  T    11
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
	size := len(kmer)
	if size == 0 || size > 32 {
		return 0, ErrKOverflow
	}

	for i := range kmer {
		switch kmer[size-1-i] {
		case 'G', 'g', 'K', 'k':
			code |= 2 << uint64(i*2)
		case 'T', 't', 'U', 'u':
			code |= 3 << uint64(i*2)
		case 'C', 'c', 'S', 's', 'B', 'b', 'Y', 'y':
			code |= 1 << uint64(i*2)
		case 'A', 'a', 'N', 'n', 'M', 'm', 'V', 'v', 'H', 'h', 'R', 'r', 'D', 'd', 'W', 'w':
			code |= 0 << uint64(i*2)
		default:
			return code, ErrIllegalBase
		}
	}
	return code, nil
}

// ErrNotConsecutiveKmers means the two Kmers are not consecutive
var ErrNotConsecutiveKmers = errors.New("unikmer: not consecutive Kmers")

// MustEncodeFromPrevious encodes from previous kmer,
// assuming the kmer and preKmer are both OK.
func MustEncodeFromPreviousKmer(kmer []byte, preKmer []byte, preCode uint64) (uint64, error) {
	preCode = preCode & ((1 << (uint(len(kmer)-1) * 2)) - 1) << 2
	switch kmer[len(kmer)-1] {
	case 'G', 'g', 'K', 'k':
		preCode += 2
	case 'T', 't', 'U', 'u':
		preCode += 3
	case 'C', 'c', 'S', 's', 'B', 'b', 'Y', 'y':
		preCode += 1
	case 'A', 'a', 'N', 'n', 'M', 'm', 'V', 'v', 'H', 'h', 'R', 'r', 'D', 'd', 'W', 'w':
		preCode += 0
	default:
		return preCode, ErrIllegalBase
	}
	return preCode, nil
}

// EncodeFromPrevious encodes from previous kmer, inspired by ntHash
func EncodeFromPreviousKmer(kmer []byte, preKmer []byte, preCode uint64) (uint64, error) {
	if len(kmer) == 0 {
		return 0, ErrKOverflow
	}
	if len(kmer) != len(preKmer) {
		return 0, ErrKMismatch
	}
	if !bytes.Equal(kmer[0:len(kmer)-1], preKmer[1:len(preKmer)]) {
		return 0, ErrNotConsecutiveKmers
	}
	return MustEncodeFromPreviousKmer(kmer, preKmer, preCode)
}

// Reverse returns code of the reversed sequence.
func Reverse(code uint64, k int) (c uint64) {
	if k <= 0 || k > 32 {
		panic(ErrKOverflow)
	}
	for i := 0; i < k; i++ {
		c |= (code & 3) << uint((k-i-1)*2)
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
		c |= (3 - (code & 3)) << uint(i*2)
		code >>= 2
	}
	return
}

// bit2base is for mapping bit to base.
var bit2base = [4]byte{'A', 'C', 'G', 'T'}

// Decode converts the code to original seq
func Decode(code uint64, k int) []byte {
	if k <= 0 || k > 32 {
		panic(ErrKOverflow)
	}
	kmer := make([]byte, k)
	for i := 0; i < k; i++ {
		kmer[k-1-i] = bit2base[code&3]
		code >>= 2
	}
	return kmer
}

// KmerCode is a struct representing a kmer in 64-bits.
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

// NewKmerCodeFromPreviousOne computes KmerCode from the previous consecutive kmer.
func NewKmerCodeFromPreviousOne(kmer []byte, preKmer []byte, preKcode KmerCode) (KmerCode, error) {
	code, err := EncodeFromPreviousKmer(kmer, preKmer, preKcode.Code)
	if err != nil {
		return KmerCode{}, err
	}
	return KmerCode{code, len(kmer)}, err
}

// NewKmerCodeMustFromPreviousOne computes KmerCode from the previous consecutive kmer,
// assuming the kmer and preKmer are both OK.
func NewKmerCodeMustFromPreviousOne(kmer []byte, preKmer []byte, preKcode KmerCode) (KmerCode, error) {
	code, err := MustEncodeFromPreviousKmer(kmer, preKmer, preKcode.Code)
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
	return kcode.Rev().Comp()
}

// Bytes returns kmer in []byte.
func (kcode KmerCode) Bytes() []byte {
	return Decode(kcode.Code, kcode.K)
}

// String returns kmer in string
func (kcode KmerCode) String() string {
	return string(Decode(kcode.Code, kcode.K))
}
