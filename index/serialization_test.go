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

package index

import (
	"bufio"
	"bytes"
	"fmt"
	"io"
	"os"
	"testing"
)

func TestIndexReadAndWrite(t *testing.T) {
	file := "test.unikidx"
	defer func() {
		err := os.Remove(file)
		if err != nil {
			t.Errorf("clean error %s", err)
		}
	}()

	k := 31
	canonical := true
	numHashes := uint8(1)
	numSigs := uint64(2)
	names := []string{"a", "b", "c", "d", "e", "f", "g", "h", "i"}
	indices := []uint32{1, 2, 3, 4, 5, 6, 7, 8, 9}
	sizes := []uint64{1, 2, 3, 4, 5, 6, 7, 8, 9}
	data := [][]byte{[]byte("aa"), []byte("bb")}
	err := write(file, k, canonical, numHashes, numSigs, names, indices, sizes, data)
	if err != nil {
		t.Errorf("write error %s", err)
	}

	reader, datas, err := read(file)
	if err != nil {
		t.Errorf("read error %s", err)
	}
	if reader.K != k {
		t.Errorf("unmatch k")
	}

	if reader.Canonical != canonical {
		t.Errorf("unmatch canonical")
	}

	if reader.NumHashes != numHashes {
		t.Errorf("unmatch NumHashes")
	}
	if reader.NumSigs != numSigs {
		t.Errorf("unmatch NumSigs")
	}
	if len(reader.Names) != len(names) {
		t.Errorf("unmatch names length")
	}
	for i, n := range names {
		if reader.Names[i] != n {
			t.Errorf("unmatch name")
		}
	}
	if len(reader.Indices) != len(indices) {
		t.Errorf("unmatch indices length")
	}
	for i, n := range indices {
		if reader.Indices[i] != n {
			t.Errorf("unmatch index")
		}
	}
	if len(reader.Sizes) != len(sizes) {
		t.Errorf("unmatch sizes length")
	}
	for i, n := range sizes {
		if reader.Sizes[i] != n {
			t.Errorf("unmatch size")
		}
	}
	if len(datas) != len(data) {
		t.Errorf("unmatch data length")
	}
	for i, d := range data {
		if bytes.Compare(d, datas[i]) != 0 {
			t.Errorf("unmatch data")
		}
	}

}

func write(file string, k int, canonical bool, numHashes uint8, numSigs uint64, names []string, indices []uint32, sizes []uint64, datas [][]byte) error {
	w, err := os.Create(file)
	if err != nil {
		return err
	}
	defer w.Close()

	outfh := bufio.NewWriter(w)
	defer outfh.Flush()

	writer, err := NewWriter(outfh, k, canonical, numHashes, numSigs, names, indices, sizes)
	if err != nil {
		return err
	}
	for _, data := range datas {
		err = writer.Write(data)
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

func read(file string) (*Reader, [][]byte, error) {
	r, err := os.Open(file)
	if err != nil {
		return nil, nil, err
	}
	defer r.Close()

	infh := bufio.NewReader(r)

	reader, err := NewReader(infh)
	if err != nil {
		return reader, nil, err
	}

	fmt.Println(reader.Header)

	datas := make([][]byte, 0, 10)
	var data []byte
	for {
		data, err = reader.Read()
		if err != nil {
			if err == io.EOF {
				break
			}
			return nil, nil, err
		}

		datas = append(datas, data)
	}

	return reader, datas, nil
}
