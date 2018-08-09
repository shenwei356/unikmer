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
	"fmt"
	"math/rand"
	"testing"
)

var randomMers [][]byte
var randomMersN = 100000

var benchMer = []byte("ACTGactgGTCAgtcaactgGTCAACTGGTCA")
var benchCode uint64
var benchKmerCode KmerCode

func init() {
	randomMers = make([][]byte, randomMersN)
	for i := 0; i < randomMersN; i++ {
		randomMers[i] = make([]byte, rand.Intn(32)+1)
		for j := range randomMers[i] {
			randomMers[i][j] = bit2base[rand.Intn(4)]
		}
	}

	// for benchmark
	var err error
	benchCode, err = Encode(benchMer)
	if err != nil {
		panic(fmt.Sprintf("init: fail to encode %s", benchMer))
	}

	benchKmerCode, err = NewKmerCode(benchMer)
	if err != nil {
		panic(fmt.Sprintf("init: fail to create KmerCode from %s", benchMer))
	}
}

// TestEncodeDecode tests encode and decode
func TestEncodeDecode(t *testing.T) {
	var kcode KmerCode
	var err error
	for _, mer := range randomMers {
		kcode, err = NewKmerCode(mer) // encode
		if err != nil {
			t.Errorf("Encode error: %s", mer)
		}

		if !bytes.Equal(mer, kcode.Bytes()) { // decode
			t.Errorf("Decode error: %s != %s ", mer, kcode.Bytes())
		}
	}
}

func TestRevComp(t *testing.T) {
	var kcode KmerCode
	for _, mer := range randomMers {
		kcode, _ = NewKmerCode(mer)

		// fmt.Printf("%s, rev:%s\n", kcode, kcode.Rev())
		if !kcode.Rev().Rev().Equal(kcode) {
			t.Errorf("Rev() error: %s, Rev(): %s", kcode, kcode.Rev())
		}
	}

	for _, mer := range randomMers {
		kcode, _ = NewKmerCode(mer)

		// fmt.Printf("%s, comp:%s\n", kcode, kcode.Comp())
		if !kcode.Comp().Comp().Equal(kcode) {
			t.Errorf("Comp() error: %s, Comp(): %s", kcode, kcode.Comp())
		}
	}
}

// BenchmarkEncode tests speed of encode()
func BenchmarkEncodeK32(b *testing.B) {
	for i := 0; i < b.N; i++ {
		Encode(benchMer)
	}
}

// BenchmarkDecode tests speed of decode
func BenchmarkDecodeK32(b *testing.B) {
	for i := 0; i < b.N; i++ {
		Decode(benchCode, len(benchMer))
	}
}

func BenchmarkRevK32(b *testing.B) {
	for i := 0; i < b.N; i++ {
		benchKmerCode.Rev()
	}
}

func BenchmarkCompK32(b *testing.B) {
	for i := 0; i < b.N; i++ {
		benchKmerCode.Comp()
	}
}

func BenchmarkRevCompK32(b *testing.B) {
	for i := 0; i < b.N; i++ {
		benchKmerCode.RevComp()
	}
}
