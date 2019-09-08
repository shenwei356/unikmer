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
)

func TestStreamVByte64(t *testing.T) {
	ntests := 10000
	tests := make([][2]uint64, ntests)
	var i int
	for ; i < ntests/4; i++ {
		tests[i] = [2]uint64{rand.Uint64(), rand.Uint64()}
	}
	for ; i < ntests/2; i++ {
		tests[i] = [2]uint64{uint64(rand.Uint32()), uint64(rand.Uint32())}
	}
	for ; i < ntests*3/4; i++ {
		tests[i] = [2]uint64{uint64(rand.Intn(65536)), uint64(rand.Intn(256))}
	}
	for ; i < ntests; i++ {
		tests[i] = [2]uint64{uint64(rand.Intn(256)), uint64(rand.Intn(256))}
	}

	for i, test := range tests {
		buf := make([]byte, 16)
		ctrl, n := PutUint64s(buf, test[0], test[1])

		result, n2 := Uint64s(ctrl, buf[0:n])
		if n2 == 0 {
			t.Errorf("#%d, wrong decoded number", i)
		}

		if result[0] != test[0] || result[1] != test[1] {
			t.Errorf("#%d, wrong decoded result: %d, %d, answer: %d, %d", i, result[0], result[1], test[0], test[1])
		}
		// fmt.Printf("%d, %d => n=%d, buf=%v\n", test[0], test[1], n, buf[0:n])
	}
}
