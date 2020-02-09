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
	"bufio"
	"fmt"
	"io"
	"os"
	"runtime"

	"github.com/shenwei356/unikmer"
	"github.com/spf13/cobra"
)

// subsetCmd represents
var subsetCmd = &cobra.Command{
	Use:   "subset",
	Short: "extract smaller k-mers from binary file",
	Long: `extract smaller k-mers from binary file

Attention:
  1. It's faster than re-counting from sequence file but in cost of losing
    few ( <= (K-k)*2 ) k-mers in the ends of sequence and its reverse complement
    sequence.

`,
	Run: func(cmd *cobra.Command, args []string) {
		opt := getOptions(cmd)
		runtime.GOMAXPROCS(opt.NumCPUs)

		var err error

		files := getFileListFromArgsAndFile(cmd, args, true, "infile-list", true)

		if len(files) > 1 {
			checkError(fmt.Errorf("no more than one file should be given"))
		}

		checkFileSuffix(extDataFile, files...)

		outFile := getFlagString(cmd, "out-prefix")
		k := getFlagPositiveInt(cmd, "kmer-len")
		if k > 32 {
			checkError(fmt.Errorf("k > 32 not supported"))
		}

		canonical := getFlagBool(cmd, "canonical")

		file := files[0]

		var infh *bufio.Reader
		var r *os.File

		infh, r, _, err = inStream(file)
		checkError(err)
		defer r.Close()

		var reader *unikmer.Reader
		reader, err = unikmer.NewReader(infh)
		checkError(err)

		if k >= reader.K {
			log.Errorf("k (%d) should be small than k size (%d) of %s", k, reader.K, file)
			return
		}

		if !isStdout(outFile) {
			outFile += extDataFile
		}
		outfh, gw, w, err := outStream(outFile, opt.Compress, opt.CompressionLevel)
		checkError(err)
		defer func() {
			outfh.Flush()
			if gw != nil {
				gw.Close()
			}
			w.Close()
		}()

		var mode uint32
		if opt.Compact {
			mode |= unikmer.UNIK_COMPACT
		}
		if canonical {
			mode |= unikmer.UNIK_CANONICAL
		}
		writer, err := unikmer.NewWriter(outfh, k, mode)
		checkError(err)

		m := make(map[uint64]struct{}, mapInitSize)

		var kcode, kcode2 unikmer.KmerCode
		var kmer []byte
		var ok bool
		var n int64
		for {
			kcode, err = reader.Read()
			if err != nil {
				if err == io.EOF {
					break
				}
				checkError(err)
			}

			kmer = kcode.Bytes()
			kmer = kmer[0:k]

			kcode2, err = unikmer.NewKmerCode(kmer)
			if err != nil {
				checkError(fmt.Errorf("fail to encode '%s': %s", kmer, err))
			}

			if canonical {
				kcode2 = kcode2.Canonical()
			}

			if _, ok = m[kcode2.Code]; !ok {
				m[kcode2.Code] = struct{}{}
				writer.WriteCode(kcode2.Code)
			}
		}

		checkError(writer.Flush())
		if opt.Verbose {
			log.Infof("%d k-mers saved to %s", n, outFile)
		}
	},
}

func init() {
	RootCmd.AddCommand(subsetCmd)

	subsetCmd.Flags().StringP("out-prefix", "o", "-", `out file prefix ("-" for stdout)`)
	subsetCmd.Flags().IntP("kmer-len", "k", 0, "k-mer length")
	subsetCmd.Flags().BoolP("canonical", "K", false, "only keep the canonical k-mers")
}
