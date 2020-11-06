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

	"github.com/cespare/xxhash"
	"github.com/shenwei356/bio/seq"
	"github.com/will-rowe/nthash"
)

// ErrInvalidS means s >= k.
var ErrInvalidS = fmt.Errorf("unikmer: invalid s-mer size")

// Sketch is a kmer sketch.
type Sketch struct {
	S        []byte
	k        int
	s        int
	l        int // 2*k-s-1
	circular bool
	hasher   *nthash.NTHi

	kMs int // k-s
	r   int // L-s

	idx int // current location, 0-based
	end int

	i, mI     int
	v, mV     uint64
	maxUint64 uint64
}

// NewSyncmerSketch returns ntHash SyncmerSketch.
func NewSyncmerSketch(S *seq.Seq, k int, s int, circular bool) (*Sketch, error) {
	if k < 1 {
		return nil, ErrInvalidK
	}
	if s >= k {
		return nil, ErrInvalidS
	}
	if len(S.Seq) < k*2-s-1 {
		return nil, ErrShortSeq
	}

	sketch := &Sketch{S: S.Seq, s: s, k: k, l: 2*k - s - 1, circular: circular}
	sketch.maxUint64 = ^uint64(0)

	var seq2 []byte
	if circular {
		seq2 = make([]byte, len(S.Seq), len(S.Seq)+k-1)
		copy(seq2, S.Seq) // do not edit original sequence
		seq2 = append(seq2, S.Seq[0:k-1]...)
		sketch.S = seq2
	} else {
		seq2 = S.Seq
	}
	sketch.end = len(seq2) - 2*k + s + 1 // len(sequence) -L
	sketch.r = 2*k - s - 1 - s           // L-s

	var err error
	sketch.hasher, err = nthash.NewHasher(&seq2, uint(k))
	if err != nil {
		return nil, err
	}

	return sketch, nil
}

// Next returns next ntHash.
func (s *Sketch) Next() (code uint64, ok bool) {

	for {
		// fmt.Println(s.idx, s.end, cap(s.S))
		if s.idx > s.end {
			return 0, false
		}

		// nthash of current k-mer
		code, ok = s.hasher.Next(true)
		if !ok {
			return code, false
		}
		// fmt.Fprintf(os.Stderr, "%d:%d-%s: %d\n", s.idx, s.idx+s.k, s.S[s.idx:s.idx+s.k], code)

		// find min s-mer
		s.mV = s.maxUint64
		for s.i = s.idx; s.i <= s.idx+s.r; s.i++ {
			s.v = xxhash.Sum64(s.S[s.i : s.i+s.s])
			// fmt.Fprintf(os.Stderr, "  %d:%d-%s: %d\n", s.i, s.i+s.s, s.S[s.i:s.i+s.s], s.v)
			if s.v < s.mV {
				s.mI, s.mV = s.i, s.v
			}
		}
		// fmt.Fprintf(os.Stderr, "  min: %d-%s\n", s.mI, s.S[s.mI:s.mI+s.s])

		// check if this k-mer is bounded syncmer
		if s.mI == s.idx || s.mI == s.idx+s.kMs { // beginning || end
			// fmt.Fprintf(os.Stderr, "  Yes: %d-%s: %d\n", s.mI, s.S[s.mI:s.mI+s.s], s.mV)
			s.idx++
			return code, true
		}

		s.idx++
	}
}

// CurrentIndex returns current  0-baesd index
func (s *Sketch) CurrentIndex() int {
	return s.idx - 1
}
