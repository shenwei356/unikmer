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

func putUvarint(buf []byte, x uint64) int {
	i := 0
	var t, j uint64

	for j = 1; j <= 7; j++ {
		t = 8 * (8 - j)
		if x > 1<<(t)-1 {
			buf[i] = byte((x & (0xff << t)) >> t)
			x = x & (1<<(8*(9-j)) - 1)
			i++
		}
	}
	buf[i] = byte(x)
	return i + 1
}

func uvarint(buf []byte, n int) uint64 {
	if n > 8 {
		panic("invalid n")
	}
	if n < 1 {
		n = len(buf)
	}
	var x uint64
	for i := n - 1; i >= 0; i-- {
		x |= uint64(buf[i]) << uint((n-1-i)*8)
	}
	return x
}
