// Copyright © 2018-2021 Wei Shen <shenwei356@gmail.com>
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

// CodeTaxid is the code-taxid pair
type CodeTaxid struct {
	Code uint64
	// _     uint32 // needed? to test
	Taxid uint32
}

// CodeTaxidSlice is a list of CodeTaxid, just for sorting
type CodeTaxidSlice []CodeTaxid

// Len return length of the slice
func (pairs CodeTaxidSlice) Len() int {
	return len(pairs)
}

// Swap swaps two elements
func (pairs CodeTaxidSlice) Swap(i, j int) {
	pairs[i], pairs[j] = pairs[j], pairs[i]
}

// Less simply compare two KmerCode
func (pairs CodeTaxidSlice) Less(i, j int) bool {
	return pairs[i].Code < pairs[j].Code
}
