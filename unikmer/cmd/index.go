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
	"io"
	"runtime"
	"strings"

	"github.com/shenwei356/unikmer"
	"github.com/shenwei356/xopen"
	"github.com/spf13/cobra"
	boom "github.com/tylertreat/BoomFilters"
)

// indexCmd represents
var indexCmd = &cobra.Command{
	Use:   "index",
	Short: "create index for binary file",
	Long: `create index for binary file

`,
	Run: func(cmd *cobra.Command, args []string) {
		opt := getOptions(cmd)
		runtime.GOMAXPROCS(opt.NumCPUs)
		files := getFileList(args)

		if len(files) == 1 && isStdin(files[0]) {
			log.Errorf("%s file needed", extDataFile)
			return
		}

		hint := getFlagPositiveInt(cmd, "esti-kmer-num")

		var infh *xopen.Reader
		var reader *unikmer.Reader
		var kcode unikmer.KmerCode
		var mer []byte
		var err error

		for _, file := range files {
			if !isStdin(file) && !strings.HasSuffix(file, extDataFile) {
				log.Errorf("input should be stdin or %s file", extDataFile)
				return
			}

			func() {
				if isStdin(file) {
					log.Warningf("no need to create index for stdin")
					return
				}

				infh, err = xopen.Ropen(file)
				checkError(err)
				defer infh.Close()

				reader, err = unikmer.NewReader(infh)
				checkError(err)

				sbf := boom.NewScalableBloomFilter(uint(hint), 0.01, 0.8)
				ibf := boom.NewInverseBloomFilter(uint(hint / 5))

				for {
					kcode, err = reader.Read()
					if err != nil {
						if err == io.EOF {
							break
						}
						checkError(err)
					}

					mer = kcode.Bytes()
					sbf.Add(mer)
					ibf.Add(mer)
				}

				// save to files
				writeIndex(reader.K, sbf, ibf, file+extSBF, file+extIBF)
			}()
		}
	},
}

func init() {
	RootCmd.AddCommand(indexCmd)

	indexCmd.Flags().IntP("esti-kmer-num", "n", 100000000, "estimated kmer num length (for initializing Bloom Filter)")

}
