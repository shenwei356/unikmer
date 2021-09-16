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

package unikmer

import (
	"bufio"
	"bytes"
	"fmt"
	"io"
	"math/rand"
	"os"
	"sort"
	"testing"

	"github.com/shenwei356/util/byteutil"
)

func genKmers(k int, num int, sorted bool) [][]byte {
	mers := make([][]byte, num)
	var j int
	for i := 0; i < num; i++ {
		mers[i] = make([]byte, k)
		for j = 0; j < k; j++ {
			mers[i][j] = bit2base[rand.Intn(4)]
		}
	}
	sort.Sort(byteutil.SliceOfByteSlice(mers))
	return mers
}

// TestWriterReader tests Writer and Writer
func TestWriter(t *testing.T) {
	var file string

	var mers, mers2 [][]byte
	var err error

	ns := []int{10001, 10001, 10001, 10000}
	for k := 1; k <= 31; k++ {
		for i, flag := range []uint32{0, UnikCompact, UnikSorted} { //, UnikSorted
			func(flag uint32) {
				mers = genKmers(k, ns[i], flag&UnikSorted > 0)

				file = fmt.Sprintf("t.k%d.unik", k)

				err = write(mers, file, flag)
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
						t.Errorf("write and read: data mismatch. %d: %d vs %d", i, mers[i], mers2[i])
					}
				}
			}(flag)
		}
	}
}

func write(mers [][]byte, file string, flag uint32) error {
	w, err := os.Create(file)
	if err != nil {
		return err
	}
	defer w.Close()

	outfh := bufio.NewWriter(w)
	defer outfh.Flush()

	writer, err := NewWriter(outfh, len(mers[0]), flag)
	if err != nil {
		return err
	}

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

	infh := bufio.NewReader(r)

	reader, err := NewReader(infh)
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
