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

import (
	"math"
)

// From https://github.com/bingmann/cobs/blob/master/cobs/util/calc_signature_size.cpp
func CalcSignatureSize(numElements uint64, numHashes int, falsePositiveRate float64) uint64 {
	ratio := float64(-numHashes) / (math.Log(1 - math.Pow(falsePositiveRate, 1/float64(numHashes))))
	return uint64(math.Ceil(float64(numElements) * ratio))
}

// get the two basic hash function values for data.
// Based on early version of https://github.com/willf/bloom/blob/master/bloom.go .
func baseHashes(key uint64) (uint32, uint32) {
	hash := hash64(key)
	return uint32(hash >> 32), uint32(hash)
}

// return locations in bitset for a key
func hashLocations(key uint64, numHashes int, numSigs uint64) []int {
	locs := make([]int, numHashes)
	a, b := baseHashes(key)
	for i := uint32(0); i < uint32(numHashes); i++ {
		locs[i] = int(uint64(a+b*i) % numSigs)
	}
	return locs
}

// return hashes for a key
func hashValues(key uint64, numHashes int) []int {
	locs := make([]int, numHashes)
	a, b := baseHashes(key)
	for i := uint32(0); i < uint32(numHashes); i++ {
		locs[i] = int(uint64(a + b*i))
	}
	return locs
}

// https://gist.github.com/badboy/6267743 .
// version with mask: https://gist.github.com/lh3/974ced188be2f90422cc .
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
