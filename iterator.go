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
	"fmt"

	"github.com/pkg/errors"
	"github.com/shenwei356/bio/seq"
	"github.com/will-rowe/nthash"
)

// ErrInvalidK means k < 1.
var ErrInvalidK = fmt.Errorf("unikmer: invalid k-mer size")

// ErrEmptySeq sequence is empty.
var ErrEmptySeq = fmt.Errorf("unikmer: empty sequence")

// ErrShortSeq means the sequence is shorter than k
var ErrShortSeq = fmt.Errorf("unikmer: sequence shorter than k")

// Iterator is a kmer code (k<=32) or hash iterator.
type Iterator struct {
	s         *seq.Seq // only used for KmerIterator
	k         int
	kUint     uint // uint(k)
	kP1       int  // k -1
	kP1Uint   uint // uint(k-1)
	canonical bool
	circular  bool

	hash bool

	finished     bool
	revcomStrand bool
	idx          int

	// for KmerIterator
	length    int
	end, e    int
	first     bool
	kmer      []byte
	codeBase  uint64
	preCode   uint64
	preCodeRC uint64
	codeRC    uint64

	// for HashIterator
	hasher *nthash.NTHi
}

// NewHashIterator returns ntHash Iterator.
func NewHashIterator(s *seq.Seq, k int, canonical bool, circular bool) (*Iterator, error) {
	if k < 1 {
		return nil, ErrInvalidK
	}
	if len(s.Seq) < k {
		return nil, ErrShortSeq
	}

	iter := &Iterator{s: s, k: k, canonical: canonical, circular: circular}
	iter.hash = true
	iter.kUint = uint(k)
	iter.kP1 = k - 1
	iter.kP1Uint = uint(k - 1)

	var err error
	var seq2 []byte
	if circular {
		seq2 = make([]byte, len(s.Seq), len(s.Seq)+k-1)
		copy(seq2, s.Seq) // do not edit original sequence
		seq2 = append(seq2, s.Seq[0:k-1]...)
	} else {
		seq2 = s.Seq
	}
	iter.hasher, err = nthash.NewHasher(&seq2, uint(k))
	if err != nil {
		return nil, err
	}

	return iter, nil
}

// NextHash returns next ntHash.
func (iter *Iterator) NextHash() (code uint64, ok bool) {
	code, ok = iter.hasher.Next(iter.canonical)
	iter.idx++
	return code, ok
}

// NewKmerIterator returns kmer code iterator
func NewKmerIterator(s *seq.Seq, k int, canonical bool, circular bool) (*Iterator, error) {
	if k < 1 {
		return nil, ErrInvalidK
	}
	if len(s.Seq) < k {
		return nil, ErrShortSeq
	}

	var s2 *seq.Seq
	if circular {
		s2 = s.Clone() // do not edit original sequence
		s2.Seq = append(s2.Seq, s.Seq[0:k-1]...)
	} else {
		s2 = s
	}

	iter := &Iterator{s: s2, k: k, canonical: canonical, circular: circular}
	iter.length = len(s2.Seq)
	iter.end = iter.length - k - 1
	iter.kUint = uint(k)
	iter.kP1 = k - 1
	iter.kP1Uint = uint(k - 1)

	iter.first = true

	return iter, nil
}

// NextKmer returns next kmer code
func (iter *Iterator) NextKmer() (code uint64, ok bool, err error) {
	if iter.finished {
		return 0, false, nil
	}

	if iter.idx == iter.end {
		if iter.canonical || iter.revcomStrand {
			iter.finished = true
			return 0, false, nil
		}
		iter.s.RevComInplace()
		iter.idx = 0
		iter.revcomStrand = true
		iter.first = true
	}

	iter.e = iter.idx + iter.k

	if iter.e > iter.length {
		if iter.canonical || iter.revcomStrand { //  end of sequence
			iter.finished = true
			return 0, false, nil
		}
	} else {
		iter.kmer = iter.s.Seq[iter.idx:iter.e]
	}

	if !iter.first {
		iter.codeBase = base2bit[iter.kmer[iter.kP1]]
		if iter.codeBase == 4 {
			err = ErrIllegalBase
		}

		// compute code from previous one
		code = iter.preCode&((1<<(iter.kP1Uint<<1))-1)<<2 | iter.codeBase

		// compute code of revcomp kmer from previous one
		iter.codeRC = (iter.codeBase^3)<<(iter.kP1Uint<<1) | (iter.preCodeRC >> 2)
	} else {
		code, err = Encode(iter.kmer)
		iter.codeRC = MustRevComp(code, iter.k)
		iter.first = false
	}
	if err != nil {
		return 0, false, errors.Wrapf(err, "encode %s", iter.kmer)
	}

	iter.preCode = code
	iter.preCodeRC = iter.codeRC
	iter.idx++

	if iter.canonical && code > iter.codeRC {
		code = iter.codeRC
	}

	return code, true, nil
}

// Next is a wrapter for NextHash and NextKmer
func (iter *Iterator) Next() (code uint64, ok bool, err error) {
	if iter.hash {
		code, ok = iter.NextHash()
		return
	}
	code, ok, err = iter.NextKmer()
	return
}

// CurrentIndex returns current  0-baesd index
func (iter *Iterator) CurrentIndex() int {
	return iter.idx - 1
}
