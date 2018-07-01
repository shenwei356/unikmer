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
	"bufio"
	"bytes"
	"compress/gzip"
	"io"
	"math/rand"
	"os"
	"testing"
)

func genKmers(k int, num int) [][]byte {
	mers := make([][]byte, num)
	var j int
	for i := 0; i < num; i++ {
		mers[i] = make([]byte, k)
		for j = 0; j < k; j++ {
			mers[i][j] = code2base[rand.Intn(4)]
		}
	}
	return mers
}

// TestWriterReader tests Writer and Writer
func TestWriter(t *testing.T) {
	var file = "t.unik.gz"

	var mers, mers2 [][]byte
	var err error

	mers = genKmers(21, 10000)

	err = write(mers, file)
	if err != nil {
		t.Error(err)
	}
	defer func() {
		err = os.Remove(file)
		if err != nil {
			t.Error(err)
		}
	}()

	mers2, err = read(file)
	if err != nil {
		t.Error(err)
	}

	if len(mers2) != len(mers) {
		t.Errorf("write and read: number err")
	}
	for i := 0; i < len(mers); i++ {
		if !bytes.Equal(mers[i], mers2[i]) {
			t.Errorf("write and read: data mismatch")
		}
	}
}

func write(mers [][]byte, file string) error {
	w, err := os.Create(file)
	if err != nil {
		return err
	}
	defer w.Close()

	bw := bufio.NewWriter(w)
	defer bw.Flush()

	gz := gzip.NewWriter(bw)
	defer gz.Close()

	writer := NewWriter(gz, len(mers[0]))
	for _, mer := range mers {
		err = writer.WriteKmer(mer)
		if err != nil {
			return err
		}
	}
	err = writer.Flush()
	if err != nil {
		return err
	}

	return nil
}

func read(file string) ([][]byte, error) {
	r, err := os.Open(file)
	if err != nil {
		return nil, err
	}
	defer r.Close()

	br := bufio.NewReader(r)

	gz, err := gzip.NewReader(br)
	if err != nil {
		return nil, err
	}
	defer gz.Close()

	reader, err := NewReader(gz)
	if err != nil {
		return nil, err
	}

	// fmt.Println(reader.Header)

	mers := make([][]byte, 0, 1000)
	var kcode KmerCode
	for {
		kcode, err = reader.Read()
		if err != nil {
			if err == io.EOF {
				break
			}
			return nil, err
		}

		mers = append(mers, kcode.Bytes())
	}

	return mers, nil
}
