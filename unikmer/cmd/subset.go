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
	Short: "extract smaller Kmers from binary file",
	Long: `extract smaller Kmers from binary file

Attention:
  - It's faster than re-counting from sequence file but in cost of losing
    few ( <= (K-k)*2 ) kmers in the ends of sequence and its reverse complement
    sequence.

`,
	Run: func(cmd *cobra.Command, args []string) {
		opt := getOptions(cmd)
		runtime.GOMAXPROCS(opt.NumCPUs)
		files := getFileList(args)

		if len(files) > 1 {
			checkError(fmt.Errorf("no more than one file should be given"))
		}

		checkFiles(files)

		outFile := getFlagString(cmd, "out-prefix")
		k := getFlagPositiveInt(cmd, "kmer-len")
		if k > 32 {
			checkError(fmt.Errorf("k > 32 not supported"))
		}

		file := files[0]

		var err error
		var infh *bufio.Reader
		var r *os.File

		infh, r, err = inStream(file)
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
		outfh, gw, w, err := outStream(outFile, opt.Compress)
		checkError(err)
		defer func() {
			outfh.Flush()
			if gw != nil {
				gw.Close()
			}
			w.Close()
		}()

		writer := unikmer.NewWriter(outfh, k)
		writer.Compact = opt.Compact

		m := make(map[uint64]struct{}, mapInitSize)

		var kcode, kcode2 unikmer.KmerCode
		var kmer []byte
		var ok bool
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
				checkError(fmt.Errorf("encoding '%s': %s", kmer, err))
			}

			if _, ok = m[kcode2.Code]; !ok {
				m[kcode2.Code] = struct{}{}
				checkError(writer.Write(kcode2))
			}
		}
	},
}

func init() {
	RootCmd.AddCommand(subsetCmd)

	subsetCmd.Flags().StringP("out-prefix", "o", "-", `out file prefix ("-" for stdout)`)
	subsetCmd.Flags().IntP("kmer-len", "k", 0, "kmer length")
}
