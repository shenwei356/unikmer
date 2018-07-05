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
	"bufio"
	"encoding/binary"
	"fmt"
	"io"
	"os"
	"sync"

	gzip "github.com/klauspost/pgzip"
	"github.com/shenwei356/unikmer"
	"github.com/shenwei356/xopen"
	boom "github.com/tylertreat/BoomFilters"
)

const extDataFile = ".unik"
const extSBF = ".sbf"
const extIBF = ".ibf"

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

func writeIndex(k int, sbf *boom.ScalableBloomFilter, ibf *boom.InverseBloomFilter, fileSBF string, fileIBF string) {
	var wg sync.WaitGroup

	wg.Add(1)
	go func() {
		defer wg.Done()
		outfhSBF, err := xopen.WopenGzip(fileSBF)
		if err != nil {
			checkError(fmt.Errorf("write %s file '%s': %s", extSBF, fileSBF, err))
		}
		defer outfhSBF.Close()

		err = writeHeader(outfhSBF, k)
		if err != nil {
			checkError(fmt.Errorf("write %s file '%s': %s", extSBF, fileSBF, err))
		}

		_, err = sbf.WriteTo(outfhSBF)
		if err != nil {
			checkError(fmt.Errorf("write %s file '%s': %s", extSBF, fileSBF, err))
		}
	}()

	wg.Add(1)
	go func() {
		defer wg.Done()
		w, err := os.Create(fileIBF)
		if err != nil {
			checkError(fmt.Errorf("write %s file '%s': %s", extIBF, fileIBF, err))
		}
		defer w.Close()

		// bw := bufio.NewWriter(w)
		// defer bw.Flush()

		outfhIBF := w
		// outfhIBF := gzip.NewWriter(w)
		// defer outfhIBF.Close()

		// outfhIBF, err := xopen.WopenGzip(fileIBF)
		// if err != nil {
		// 	checkError(fmt.Errorf("write %s file '%s': %s", extIBF, fileIBF, err))
		// }
		// defer outfhIBF.Close()

		err = writeHeader(outfhIBF, k)
		if err != nil {
			checkError(fmt.Errorf("write %s file '%s': %s", extIBF, fileIBF, err))
		}

		_, err = ibf.WriteTo(outfhIBF)
		if err != nil {
			checkError(fmt.Errorf("write %s file '%s': %s", extIBF, fileIBF, err))
		}
	}()

	wg.Wait()
}

func readIndex(fileSBF string, fileIBF string) (sbf *boom.ScalableBloomFilter, ibf *boom.InverseBloomFilter, k int) {
	var headerSBF, headerIBF unikmer.Header

	var wg sync.WaitGroup

	wg.Add(1)
	go func() {
		defer wg.Done()

		r, err := os.Open(fileSBF)
		if err != nil {
			checkError(fmt.Errorf("read %s file '%s': %s", extSBF, fileSBF, err))
		}
		defer r.Close()

		br := bufio.NewReader(r)

		var infh *gzip.Reader
		infh, err = gzip.NewReader(br)
		if err != nil {
			checkError(fmt.Errorf("read %s file '%s': %s", extSBF, fileSBF, err))
		}
		defer infh.Close()

		headerSBF, err = readHeader(infh)
		if err != nil {
			checkError(fmt.Errorf("read %s file '%s': %s", extSBF, fileSBF, err))
		}

		sbf = &boom.ScalableBloomFilter{}
		_, err = sbf.ReadFrom(infh)
		if err != nil {
			checkError(fmt.Errorf("read %s file '%s': %s", extSBF, fileSBF, err))
		}
	}()

	wg.Add(1)
	go func() {
		defer wg.Done()

		r, err := os.Open(fileIBF)
		if err != nil {
			checkError(fmt.Errorf("read %s file '%s': %s", extIBF, fileIBF, err))
		}
		defer r.Close()

		// br := bufio.NewReader(r)

		// infh, err := gzip.NewReader(br)
		// if err != nil {
		// 	checkError(fmt.Errorf("read %s file '%s': %s", extIBF, fileIBF, err))
		// }
		// defer infh.Close()

		infh := r

		headerIBF, err = readHeader(infh)
		if err != nil {
			checkError(fmt.Errorf("read %s file '%s': %s", extIBF, fileIBF, err))
		}

		ibf = boom.NewInverseBloomFilter(10000)
		_, err = ibf.ReadFrom(infh)
		if err != nil {
			checkError(fmt.Errorf("read %s file '%s': %s", extIBF, fileIBF, err))
		}

	}()

	wg.Wait()

	if headerIBF.Version != headerSBF.Version {
		checkError(fmt.Errorf("version mismatch: %s (%s) != %s (%s)", headerSBF.Version, extSBF, headerIBF.Version, extIBF))
	}
	if headerIBF.K != headerSBF.K {
		checkError(fmt.Errorf("k size mismatch: %d (%s) != %d (%s)", headerSBF.K, extSBF, headerIBF.K, extIBF))
	}

	k = headerIBF.K

	return
}
