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
	"errors"
	"fmt"
	"io"
	"os"

	gzip "github.com/klauspost/pgzip"
)

func outStream(file string, gzipped bool) (*bufio.Writer, io.WriteCloser, *os.File, error) {
	var err error
	var w *os.File
	if file == "-" {
		w = os.Stdout
	} else {
		w, err = os.Create(file)
		if err != nil {
			return nil, nil, nil, fmt.Errorf("fail to write %s: %s", file, err)
		}
	}

	if gzipped {
		gw := gzip.NewWriter(w)
		return bufio.NewWriterSize(gw, os.Getpagesize()), gw, w, nil
	}
	return bufio.NewWriterSize(w, os.Getpagesize()), nil, w, nil
}

func inStream(file string) (*bufio.Reader, *os.File, error) {
	var err error
	var r *os.File
	if file == "-" {
		if !detectStdin() {
			return nil, nil, errors.New("stdin not detected")
		}
		r = os.Stdin
	} else {
		r, err = os.Open(file)
		if err != nil {
			return nil, nil, fmt.Errorf("fail to read %s: %s", file, err)
		}
	}

	br := bufio.NewReaderSize(r, os.Getpagesize())

	if gzipped, err := isGzip(br); err != nil {
		return nil, nil, fmt.Errorf("fail to check is file (%s) gzipped: %s", file, err)
	} else if gzipped {
		gr, err := gzip.NewReader(br)
		if err != nil {
			return nil, r, fmt.Errorf("fail to create gzip reader for %s: %s", file, err)
		}
		br = bufio.NewReaderSize(gr, os.Getpagesize())
	}

	return br, r, nil
}

func isGzip(b *bufio.Reader) (bool, error) {
	return checkBytes(b, []byte{0x1f, 0x8b})
}

func checkBytes(b *bufio.Reader, buf []byte) (bool, error) {
	m, err := b.Peek(len(buf))
	if err != nil {
		return false, fmt.Errorf("no content")
	}
	for i := range buf {
		if m[i] != buf[i] {
			return false, nil
		}
	}
	return true, nil
}

func detectStdin() bool {
	// http://stackoverflow.com/a/26567513
	stat, err := os.Stdin.Stat()
	if err != nil {
		return false
	}
	return (stat.Mode() & os.ModeCharDevice) == 0
}