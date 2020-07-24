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

package cmd

// https://gist.github.com/badboy/6267743
func hash64(key uint64) uint64 {
	key = (^key) + (key << 21) // key = (key << 21) - key - 1
	key = key ^ (key >> 24)
	key = (key + (key << 3)) + (key << 8) // key * 265
	key = key ^ (key >> 14)
	key = (key + (key << 2)) + (key << 4) // key * 21
	key = key ^ (key >> 28)
	key = key + (key << 31)
	return key
}

// https://naml.us/post/inverse-of-a-hash-function/
func ihash64(key uint64) uint64 {
	var tmp uint64

	// Invert key = key + (key << 31)
	tmp = key - (key << 31)
	key = key - (tmp << 31)

	// Invert key = key ^ (key >> 28)
	tmp = key ^ key>>28
	key = key ^ tmp>>28

	// Invert key *= 21
	key *= 14933078535860113213

	// Invert key = key ^ (key >> 14)
	tmp = key ^ key>>14
	tmp = key ^ tmp>>14
	tmp = key ^ tmp>>14
	key = key ^ tmp>>14

	// Invert key *= 265
	key *= 15244667743933553977

	// Invert key = key ^ (key >> 24)
	tmp = key ^ key>>24
	key = key ^ tmp>>24

	// Invert key = (^key) + (key << 21)
	tmp = ^key
	tmp = ^(key - (tmp << 21))
	tmp = ^(key - (tmp << 21))
	key = ^(key - (tmp << 21))

	return key
}
