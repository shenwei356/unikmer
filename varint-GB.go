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

var offsets = []uint8{56, 48, 40, 32, 24, 16, 8, 0}

// PutUint64s endcodes two uint64s into 2-16 bytes, and returns control byte
// and encoded byte length.
func PutUint64s(buf []byte, v1, v2 uint64) (ctrl byte, n int) {
	blen := byteLength(v1)
	ctrl |= byte(blen - 1)
	for _, offset := range offsets[8-blen:] {
		buf[n] = byte((v1 >> offset) & 0xff)
		n++
	}

	ctrl <<= 3
	blen = byteLength(v2)
	ctrl |= byte(blen - 1)
	for _, offset := range offsets[8-blen:] {
		buf[n] = byte((v2 >> offset) & 0xff)
		n++
	}
	return
}

// Uint64s decode from encoded bytes
func Uint64s(ctrl byte, buf []byte) (values [2]uint64, n int) {
	blens := ctrlByte2ByteLengths[ctrl]
	if len(buf) < int(blens[0]+blens[1]) {
		return values, 0
	}
	for i := 0; i < 2; i++ {
		for j := uint8(0); j < blens[i]; j++ {
			values[i] <<= 8
			values[i] |= uint64(buf[n])
			n++
		}
	}

	return
}

func byteLength(n uint64) uint8 {
	if n < 256 {
		return 1
	}
	if n < 65536 {
		return 2
	}
	if n < 16777216 {
		return 3
	}
	if n < 4294967296 {
		return 4
	}
	if n < 1099511627776 {
		return 5
	}
	if n < 281474976710656 {
		return 6
	}
	if n < 72057594037927936 {
		return 7
	}
	return 8
}

var ctrlByte2ByteLengths = [64][2]uint8{
	[2]uint8{1, 1}, // 0, 0b000000
	[2]uint8{1, 2},
	[2]uint8{1, 3},
	[2]uint8{1, 4},
	[2]uint8{1, 5},
	[2]uint8{1, 6},
	[2]uint8{1, 7},
	[2]uint8{1, 8},
	[2]uint8{2, 1}, // 8, 0b001000
	[2]uint8{2, 2},
	[2]uint8{2, 3},
	[2]uint8{2, 4},
	[2]uint8{2, 5},
	[2]uint8{2, 6},
	[2]uint8{2, 7},
	[2]uint8{2, 8},
	[2]uint8{3, 1}, // 16, 0b010000
	[2]uint8{3, 2},
	[2]uint8{3, 3},
	[2]uint8{3, 4},
	[2]uint8{3, 5},
	[2]uint8{3, 6},
	[2]uint8{3, 7},
	[2]uint8{3, 8},
	[2]uint8{4, 1}, // 24, 0b011000
	[2]uint8{4, 2},
	[2]uint8{4, 3},
	[2]uint8{4, 4},
	[2]uint8{4, 5},
	[2]uint8{4, 6},
	[2]uint8{4, 7},
	[2]uint8{4, 8},
	[2]uint8{5, 1}, // 32, 0b100000
	[2]uint8{5, 2},
	[2]uint8{5, 3},
	[2]uint8{5, 4},
	[2]uint8{5, 5},
	[2]uint8{5, 6},
	[2]uint8{5, 7},
	[2]uint8{5, 8},
	[2]uint8{6, 1}, // 40, 0b101000
	[2]uint8{6, 2},
	[2]uint8{6, 3},
	[2]uint8{6, 4},
	[2]uint8{6, 5},
	[2]uint8{6, 6},
	[2]uint8{6, 7},
	[2]uint8{6, 8},
	[2]uint8{7, 1}, // 48, 0b110000
	[2]uint8{7, 2},
	[2]uint8{7, 3},
	[2]uint8{7, 4},
	[2]uint8{7, 5},
	[2]uint8{7, 6},
	[2]uint8{7, 7},
	[2]uint8{7, 8},
	[2]uint8{8, 1}, // 56, 0b111000
	[2]uint8{8, 2},
	[2]uint8{8, 3},
	[2]uint8{8, 4},
	[2]uint8{8, 5},
	[2]uint8{8, 6},
	[2]uint8{8, 7},
	[2]uint8{8, 8},
}
