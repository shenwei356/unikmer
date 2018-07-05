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

package cmd

import (
	"encoding/binary"
	"fmt"
	"io"

	"github.com/shenwei356/unikmer"
)

const extDataFile = ".unik"

var be = binary.BigEndian

func readHeader(r io.Reader) (unikmer.Header, error) {
	// check Magic number
	var m [8]byte
	var err error
	err = binary.Read(r, be, &m)
	if err != nil {
		return unikmer.Header{}, err
	}
	same := true
	for i := 0; i < 8; i++ {
		if unikmer.Magic[i] != m[i] {
			same = false
			break
		}
	}
	if !same {
		return unikmer.Header{}, unikmer.ErrInvalidFileFormat
	}

	// read metadata
	var meta [3]int64
	err = binary.Read(r, be, &meta)
	if err != nil {
		return unikmer.Header{}, err
	}
	return unikmer.Header{Version: fmt.Sprintf("%d.%d", meta[0], meta[1]),
		K: int(meta[2])}, nil
}

func writeHeader(w io.Writer, k int) error {
	err := binary.Write(w, be, unikmer.Magic)
	if err != nil {
		return err
	}
	err = binary.Write(w, be, [3]int64{unikmer.MainVersion, unikmer.MinorVersion, int64(k)})
	if err != nil {
		return err
	}
	return nil
}
