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
	"math/rand"
	"testing"

	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/util/bytesize"
)

func TestSyncmer(t *testing.T) {
	_s := "GGCAAGTTCGTCATCGATC"
	// _s := "GGCAAGTTC"
	sequence, err := seq.NewSeq(seq.DNA, []byte(_s))
	if err != nil {
		t.Errorf("fail to create sequence: %s", _s)
	}
	k := 5
	s := 2

	sketch, err := NewSyncmerSketch(sequence, k, s, false)
	if err != nil {
		t.Errorf("fail to create syncmer sketch")
	}
	var code uint64
	var ok bool
	var idx int
	codes := make([]uint64, 0, 1024)
	for {
		code, ok = sketch.NextSyncmer()
		// fmt.Println(sketch.Index(), code, ok)
		if !ok {
			break
		}

		idx = sketch.Index()

		_syncmerIdx = idx
		_syncmer = code

		codes = append(codes, code)
		// fmt.Printf("syncmer: %d-%s, %d\n", idx, _s[idx:idx+k], code)
	}
	// if len(codes) == 5 &&
	// 	codes[0] == 7385093395039290540 &&
	// 	codes[1] == 1099502864234245338 {
	// } else {
	// 	t.Errorf("syncmer error")
	// }
}

func TestMinimizer(t *testing.T) {
	_s := "GGCAAGTTCGTCA"
	// _s := "GGCAAGTTC"
	sequence, err := seq.NewSeq(seq.DNA, []byte(_s))
	if err != nil {
		t.Errorf("fail to create sequence: %s", _s)
	}
	k := 5
	w := 3

	sketch, err := NewMinimizerSketch(sequence, k, w, false)
	if err != nil {
		t.Errorf("fail to create minizimer sketch")
	}

	var code uint64
	var ok bool
	var idx int
	codes := make([]uint64, 0, 1024)
	for {
		code, ok = sketch.NextMinimizer()
		if !ok {
			break
		}

		idx = sketch.Index()

		_syncmerIdx = idx
		_syncmer = code

		codes = append(codes, code)
		// fmt.Printf("minizimer: %d-%s, %d\n", idx, _s[idx:idx+k], code)
	}

	if len(codes) == 4 &&
		codes[0] == 2645801399420473919 &&
		codes[1] == 1099502864234245338 &&
		codes[2] == 6763474888237448943 &&
		codes[3] == 2737971715116251183 {
	} else {
		t.Errorf("syncmer error")
	}
}

var _syncmer, _code uint64
var _syncmerIdx int

var benchSeqs []*seq.Seq

func init() {
	rand.Seed(11)

	sizes := []int{1 << 10} //, 1 << 20, 10 << 20}
	benchSeqs = make([]*seq.Seq, len(sizes))
	var err error
	for i, size := range sizes {
		sequence := make([]byte, size)

		// fmt.Printf("generating pseudo DNA with lenght of %s ...\n", bytesize.ByteSize(size))
		for j := 0; j < size; j++ {
			sequence[j] = bit2base[rand.Intn(4)]
		}
		benchSeqs[i], err = seq.NewSeq(seq.DNA, sequence)
		if err != nil {
			panic("should not happen")
		}
		// fmt.Println(benchSeqs[i])
	}
	// fmt.Printf("%d DNA sequences generated\n", len(sizes))
}

func BenchmarkHashIterator(b *testing.B) {
	for i := range benchSeqs {
		size := len(benchSeqs[i].Seq)
		b.Run(bytesize.ByteSize(size).String(), func(b *testing.B) {
			var code uint64
			var ok bool

			for j := 0; j < b.N; j++ {
				iter, err := NewHashIterator(benchSeqs[i], 31, true, false)
				if err != nil {
					b.Errorf("fail to create hash iterator. seq length: %d", size)
				}

				for {
					code, ok = iter.NextHash()
					if !ok {
						break
					}

					_code = code
				}
			}
		})
	}
}

func BenchmarkKmerIterator(b *testing.B) {
	for i := range benchSeqs {
		size := len(benchSeqs[i].Seq)
		b.Run(bytesize.ByteSize(size).String(), func(b *testing.B) {
			var code uint64
			var ok bool

			for j := 0; j < b.N; j++ {
				iter, err := NewKmerIterator(benchSeqs[i], 31, true, false)
				if err != nil {
					b.Errorf("fail to create hash iterator. seq length: %d", size)
				}
				for {
					code, ok, err = iter.NextKmer()
					if err != nil {
						b.Errorf("fail to get kmer code: %d-%s", iter.Index(),
							benchSeqs[i].Seq[iter.Index():iter.Index()+31])
					}

					if !ok {
						break
					}

					_code = code
				}
			}
		})
	}
}

func BenchmarkMinimizerIterator(b *testing.B) {
	for i := range benchSeqs {
		size := len(benchSeqs[i].Seq)
		b.Run(bytesize.ByteSize(size).String(), func(b *testing.B) {
			var code uint64
			var ok bool
			// var n int

			for j := 0; j < b.N; j++ {
				iter, err := NewMinimizerSketch(benchSeqs[i], 31, 15, false)
				if err != nil {
					b.Errorf("fail to create minizimer sketch. seq length: %d", size)
				}

				// n = 0
				for {
					code, ok = iter.NextMinimizer()
					if !ok {
						break
					}

					// fmt.Printf("minizimer: %d-%d\n", iter.Index(), code)

					_code = code
					// n++
				}

			}
			// fmt.Printf("minizimer for %s DNA, c=%.6f\n", bytesize.ByteSize(size).String(), float64(size)/float64(n))
		})
	}
}

// go test -v -test.bench=BenchmarkSyncmerIterator -cpuprofile profile.out -test.run=damnit
// go tool pprof -http=:8080 profile.out
func BenchmarkSyncmerIterator(b *testing.B) {
	for i := range benchSeqs {
		size := len(benchSeqs[i].Seq)
		b.Run(bytesize.ByteSize(size).String(), func(b *testing.B) {
			var code uint64
			var ok bool
			// var n int

			for j := 0; j < b.N; j++ {
				iter, err := NewSyncmerSketch(benchSeqs[i], 31, 16, false)
				if err != nil {
					b.Errorf("fail to create syncmer sketch. seq length: %d", size)
				}

				// n = 0
				for {
					code, ok = iter.NextSyncmer()
					if !ok {
						break
					}

					// fmt.Printf("syncmer: %d-%d\n", iter.Index(), code)

					_code = code
					// n++
				}

			}
			// fmt.Printf("syncmer for %s DNA, c=%.6f\n", bytesize.ByteSize(size).String(), float64(size)/float64(n))
		})
	}
}
