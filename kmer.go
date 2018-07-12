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
	"errors"
)

// ErrIllegalBase means that base beyond "ACGTU" was detected
var ErrIllegalBase = errors.New("unikmer: illegal base")

// ErrKOverflow means K > 32
var ErrKOverflow = errors.New("unikmer: K (1-32) overflow")

// Encode converts byte slice to bits.
//
//    M       AC
//    V       ACG
//    H       ACT
//    R       AG
//    D       AGT
//    W       AT
//    S       CG
//    B       CGT
//    Y       CT
//    K       GT
//
//
func Encode(mer []byte) (code uint64, err error) {
	size := len(mer)
	if size == 0 || size > 32 {
		return 0, ErrKOverflow
	}
	for i := range mer {
		switch mer[size-1-i] {
		case 'A', 'a', 'N', 'n', 'M', 'm', 'V', 'v', 'H', 'h', 'R', 'r', 'D', 'd', 'W', 'w':
			code += uint64(0 << uint(i*2))
		case 'C', 'c', 'S', 's', 'B', 'b', 'Y', 'y':
			code += uint64(1 << uint(i*2))
		case 'G', 'g', 'K', 'k':
			code += uint64(2 << uint(i*2))
		case 'T', 't', 'U', 'u':
			code += uint64(3 << uint(i*2))
		default:
			return code, ErrIllegalBase
		}
	}
	return code, nil
}

// Reverse returns code of reversed sequence
func Reverse(code uint64, k int) (c uint64) {
	if k <= 0 || k > 32 {
		panic(ErrKOverflow)
	}
	for i := 0; i < k; i++ {
		c += (code % 4) << uint((k-i-1)*2)
		code >>= 2
	}
	return
}

// Complement return code of complement sequence
func Complement(code uint64, k int) (c uint64) {
	if k <= 0 || k > 32 {
		panic(ErrKOverflow)
	}
	for i := 0; i < k; i++ {
		c += (3 - (code % 4)) << uint(i*2)
		code >>= 2
	}
	return
}

// code2base is for mapping code to base
var code2base = [4]byte{'A', 'C', 'G', 'T'}

// Decode converts the bits to origional seq
func Decode(code uint64, k int) []byte {
	if k <= 0 || k > 32 {
		panic(ErrKOverflow)
	}
	mer := make([]byte, k)
	for i := 0; i < k; i++ {
		mer[k-1-i] = code2base[code&3]
		code >>= 2
	}
	return mer
}

// KmerCode is a struct representing a kmer in 64-bits
type KmerCode struct {
	Code uint64
	K    int
}

// NewKmerCode returns a new KmerCode from byte slice
func NewKmerCode(mer []byte) (KmerCode, error) {
	code, err := Encode(mer)
	if err != nil {
		return KmerCode{}, err
	}
	return KmerCode{code, len(mer)}, err
}

// Equal checks wether two KmerCodes are the same
func (kcode KmerCode) Equal(kcode2 KmerCode) bool {
	return kcode.K == kcode2.K && kcode.Code == kcode2.Code
}

// Rev returns KmerCode of the reverse sequence
func (kcode KmerCode) Rev() KmerCode {
	return KmerCode{Reverse(kcode.Code, kcode.K), kcode.K}
}

// Comp returns KmerCode of the complement sequence
func (kcode KmerCode) Comp() KmerCode {
	return KmerCode{Complement(kcode.Code, kcode.K), kcode.K}
}

// RevComp returns KmerCode of the reverse complement sequence
func (kcode KmerCode) RevComp() KmerCode {
	return kcode.Rev().Comp()
}

// Bytes returns kmer in []byte
func (kcode KmerCode) Bytes() []byte {
	return Decode(kcode.Code, kcode.K)
}

// String returns kmer in string
func (kcode KmerCode) String() string {
	return string(Decode(kcode.Code, kcode.K))
}
