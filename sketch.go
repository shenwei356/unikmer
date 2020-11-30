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
	"sort"

	"github.com/cespare/xxhash"
	"github.com/shenwei356/bio/seq"
	"github.com/will-rowe/nthash"
)

// ErrInvalidS means s >= k.
var ErrInvalidS = fmt.Errorf("unikmer: invalid s-mer size")

// ErrInvalidW means w < 2 or w > (1<<32)-1
var ErrInvalidW = fmt.Errorf("unikmer: invalid minimimzer window")

// Sketch is a k-mer sketch iterator
type Sketch struct {
	S        []byte
	k        int
	s        int
	circular bool
	hasher   *nthash.NTHi

	kMs int // k-s
	r   int // L-s

	idx int // current location, 0-based
	end int

	i, mI int
	v, mV uint64

	buf     []idxValue
	i2v     idxValue
	flag    bool
	t, b, e int

	// ------ for minimizer -----
	skip      bool
	minimizer bool
	w         int
	l         int // k+w-1
	preMI     uint64
}

// NewMinimizerSketch returns a SyncmerSketch Iterator.
// It returns the minHashes in all windows of w (w>=1) bp.
func NewMinimizerSketch(S *seq.Seq, k int, w int, circular bool) (*Sketch, error) {
	if k < 1 {
		return nil, ErrInvalidK
	}
	if w < 1 || w > (1<<32)-1 {
		return nil, ErrInvalidW
	}
	if len(S.Seq) < k+w-1 {
		return nil, ErrShortSeq
	}

	sketch := &Sketch{S: S.Seq, w: w, k: k, l: k + w - 1, circular: circular}
	sketch.minimizer = true
	sketch.skip = w == 1

	var seq2 []byte
	if circular {
		seq2 = make([]byte, len(S.Seq), len(S.Seq)+k-1)
		copy(seq2, S.Seq) // do not edit original sequence
		seq2 = append(seq2, S.Seq[0:k-1]...)
		sketch.S = seq2
	} else {
		seq2 = S.Seq
	}
	sketch.end = len(seq2) - 1
	sketch.r = w - 1 // L-k

	var err error
	sketch.hasher, err = nthash.NewHasher(&seq2, uint(k))
	if err != nil {
		return nil, err
	}

	sketch.buf = make([]idxValue, 0, sketch.r+1)

	return sketch, nil
}

// NewSyncmerSketch returns a SyncmerSketch Iterator.
// 1<=s<=k.
func NewSyncmerSketch(S *seq.Seq, k int, s int, circular bool) (*Sketch, error) {
	if k < 1 {
		return nil, ErrInvalidK
	}
	if s > k || s == 0 {
		return nil, ErrInvalidS
	}
	if len(S.Seq) < k*2-s-1 {
		return nil, ErrShortSeq
	}

	sketch := &Sketch{S: S.Seq, s: s, k: k, circular: circular}
	// sketch.maxUint64 = ^uint64(0)
	sketch.skip = s == k

	var seq2 []byte
	if circular {
		seq2 = make([]byte, len(S.Seq), len(S.Seq)+k-1)
		copy(seq2, S.Seq) // do not edit original sequence
		seq2 = append(seq2, S.Seq[0:k-1]...)
		sketch.S = seq2
	} else {
		seq2 = S.Seq
	}
	sketch.end = len(seq2) - 2*k + s + 1 // len(sequence) - L (2*k - s - 1)
	sketch.r = 2*k - s - 1 - s           // L-s

	var err error
	sketch.hasher, err = nthash.NewHasher(&seq2, uint(k))
	if err != nil {
		return nil, err
	}

	sketch.buf = make([]idxValue, 0, sketch.r+1)

	return sketch, nil
}

// NextMinimizer returns next minimizer.
func (s *Sketch) NextMinimizer() (code uint64, ok bool) {
	for {
		if s.idx > s.end {
			return 0, false
		}

		// nthash of current k-mer
		code, ok = s.hasher.Next(true)
		if !ok {
			return code, false
		}

		if s.skip {
			s.mI = s.idx
			s.idx++
			return code, true
		}
		// find min k-mer
		if s.idx > s.r {
			// remove k-mer not in this window.
			// have to check position/index one by one
			for s.i, s.i2v = range s.buf {
				if s.i2v.idx == s.idx-s.w {
					if s.i < s.r {
						copy(s.buf[s.i:s.r], s.buf[s.i+1:])
					} // happen to be at the end
					s.buf = s.buf[:s.r]
					break
				}
			}

			// add new k-mer
			s.flag = false
			// using binary search, faster han linnear search
			s.b, s.e = 0, s.r-1
			for {
				s.t = s.b + (s.e-s.b)/2
				if code < s.buf[s.t].val {
					s.e = s.t - 1 // end search here
					if s.e <= s.b {
						s.flag = true
						s.i = s.b
						break
					}
				} else {
					s.b = s.t + 1 // start here
					if s.b >= s.r {
						s.flag = false
						break
					}
					if s.b >= s.e {
						s.flag = true
						s.i = s.e // right here
						break
					}
				}
			}
			if !s.flag { // it's the biggest one, append to the end
				s.buf = append(s.buf, idxValue{s.idx, code})
			} else {
				if code >= s.buf[s.i].val { // have to check again
					s.i++
				}
				s.buf = append(s.buf, blankI2V)     // append one element
				copy(s.buf[s.i+1:], s.buf[s.i:s.r]) // move right
				s.buf[s.i] = idxValue{s.idx, code}
			}

			s.i2v = s.buf[0]
			if s.i2v.val == s.preMI { // deduplicate
				s.idx++
				continue
			}

			s.mI, s.mV = s.i2v.idx, s.i2v.val
			s.preMI = s.mV

			s.idx++
			return s.i2v.val, true
		}

		if s.idx == s.r { // position w
			s.buf = append(s.buf, idxValue{idx: s.idx, val: code})
			sort.Sort(idxValues(s.buf)) // sort
		} else { // front of w
			s.buf = append(s.buf, idxValue{idx: s.idx, val: code})
		}

		s.idx++
	}
}

// NextSyncmer returns next syncmer.
func (s *Sketch) NextSyncmer() (code uint64, ok bool) {
	for {
		if s.idx > s.end {
			return 0, false
		}

		// nthash of current k-mer
		code, ok = s.hasher.Next(true)
		if !ok {
			return code, false
		}

		if s.skip {
			s.idx++
			return code, true
		}

		// find min s-mer
		if s.idx == 0 {
			for s.i = s.idx; s.i <= s.idx+s.r; s.i++ {
				s.v = xxhash.Sum64(s.S[s.i : s.i+s.s])
				s.buf = append(s.buf, idxValue{idx: s.i, val: s.v})
			}
			sort.Sort(idxValues(s.buf))
		} else {
			// remove s-mer not in this window.
			// have to check position/index one by one
			for s.i, s.i2v = range s.buf {
				if s.i2v.idx == s.idx-1 {
					if s.i < s.r {
						copy(s.buf[s.i:s.r], s.buf[s.i+1:])
					} // happen to be at the end
					s.buf = s.buf[:s.r]
					break
				}
			}

			// add new s-mer
			s.v = xxhash.Sum64(s.S[s.idx+s.r : s.idx+s.r+s.s])
			s.flag = false
			// using binary search, faster han linnear search
			s.b, s.e = 0, s.r-1
			for {
				s.t = s.b + (s.e-s.b)/2
				if s.v < s.buf[s.t].val {
					s.e = s.t - 1 // end search here
					if s.e <= s.b {
						s.flag = true
						s.i = s.b
						break
					}
				} else {
					s.b = s.t + 1 // start here
					if s.b >= s.r {
						s.flag = false
						break
					}
					if s.b >= s.e {
						s.flag = true
						s.i = s.e // right here
						break
					}
				}
			}
			if !s.flag { // it's the biggest one, append to the end
				s.buf = append(s.buf, idxValue{s.idx + s.r, s.v})
			} else {
				if s.v >= s.buf[s.i].val { // have to check again
					s.i++
				}
				s.buf = append(s.buf, blankI2V)     // append one element
				copy(s.buf[s.i+1:], s.buf[s.i:s.r]) // move right
				s.buf[s.i] = idxValue{s.idx + s.r, s.v}
			}

		}

		s.i2v = s.buf[0]
		s.mI, s.mV = s.i2v.idx, s.i2v.val

		// check if this k-mer is bounded syncmer
		if s.mI == s.idx || s.mI == s.idx+s.kMs { // beginning || end
			s.idx++
			return code, true
		}

		s.idx++
	}
}

// Next returns next sketch
func (s *Sketch) Next() (uint64, bool) {
	if s.minimizer {
		return s.NextMinimizer()
	}
	return s.NextSyncmer()
}

// Index returns current  0-baesd index
func (s *Sketch) Index() int {
	if s.minimizer {
		return s.mI
	}
	return s.idx - 1
}

// for sorting s-mer
type idxValue struct {
	idx int    // index
	val uint64 // hash
}

var blankI2V = idxValue{0, 0}

type idxValues []idxValue

func (l idxValues) Len() int               { return len(l) }
func (l idxValues) Less(i int, j int) bool { return l[i].val < l[j].val }
func (l idxValues) Swap(i int, j int)      { l[i], l[j] = l[j], l[i] }
