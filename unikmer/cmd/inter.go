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
	"fmt"
	"io"
	"runtime"
	"strings"

	"github.com/shenwei356/unikmer"
	"github.com/shenwei356/xopen"
	"github.com/spf13/cobra"
	boom "github.com/tylertreat/BoomFilters"
)

// interCmd represents
var interCmd = &cobra.Command{
	Use:   "inter",
	Short: "intersection of multiple binary files",
	Long: `intersection of multiple binary files

`,
	Run: func(cmd *cobra.Command, args []string) {
		opt := getOptions(cmd)
		runtime.GOMAXPROCS(opt.NumCPUs)
		files := getFileList(args)

		if len(files) <= 1 {
			checkError(fmt.Errorf("at least one file should be given"))
		}

		outFile := getFlagString(cmd, "out-prefix")
		hint := getFlagPositiveInt(cmd, "esti-kmer-num")

		var err error

		var ibf, ibf2 *boom.InverseBloomFilter
		// var sbf, sbf2 *boom.ScalableBloomFilter

		// sbf = boom.NewScalableBloomFilter(uint(hint), 0.01, 0.8)
		ibf = boom.NewInverseBloomFilter(uint(hint / 5))

		var infh *xopen.Reader
		var reader *unikmer.Reader
		var kcode unikmer.KmerCode
		var k int = -1
		var mer []byte
		var firstFile bool = true
		var lastFile bool
		var hasInter bool
		for i, file := range files {
			if !isStdin(file) && !strings.HasSuffix(file, extDataFile) {
				log.Errorf("input should be stdin or %s file", extDataFile)
				return
			}

			if opt.Verbose {
				log.Infof("read kmers from %s", file)
			}

			if !firstFile {
				// sbf2 = boom.NewScalableBloomFilter(uint(hint), 0.01, 0.8)
				ibf2 = boom.NewInverseBloomFilter(uint(hint / 5))
			}
			if i == len(files)-1 {
				lastFile = true
			}

			infh, err = xopen.Ropen(file)
			checkError(err)
			defer infh.Close()

			reader, err = unikmer.NewReader(infh)
			checkError(err)

			if k == -1 {
				k = reader.K
			} else if k != reader.K {
				checkError(fmt.Errorf("K (%d) of binary file '%s' not equal to previous K (%d)", reader.K, file, k))
			}
			for {
				kcode, err = reader.Read()
				if err != nil {
					if err == io.EOF {
						break
					}
					checkError(err)
				}

				mer = kcode.Bytes()

				if firstFile {
					ibf.Add(mer)
					continue
				}

				if ibf.Test(mer) { // in previous file
					ibf2.Add(mer) // add intersection into new inverse bloom filter
					if lastFile {
						hasInter = true
					}
				} else {

				}
			}

			if firstFile {
				firstFile = false
			} else {
				ibf = ibf2
				// sbf = sbf2
			}
		}

		if !hasInter {
			if opt.Verbose {
				log.Infof("no intersection found")
			}
			return
		}
		// output result from last file

		file := files[len(files)-1]
		if opt.Verbose {
			log.Infof("export kmers from last file: %s", file)
		}

		outfh, err := xopen.Wopen(outFile + extDataFile)
		checkError(err)
		defer outfh.Close()

		writer := unikmer.NewWriter(outfh, k)

		infh, err = xopen.Ropen(file)
		checkError(err)
		defer infh.Close()

		reader, err = unikmer.NewReader(infh)
		checkError(err)

		var n int64
		for {
			kcode, err = reader.Read()
			if err != nil {
				if err == io.EOF {
					break
				}
				checkError(err)
			}

			mer = kcode.Bytes()

			if ibf2.Test(mer) {
				n++
				writer.Write(kcode)
			}
		}
		if opt.Verbose {
			log.Infof("intersection of %d kmers saved", n)
		}
	},
}

func init() {
	RootCmd.AddCommand(interCmd)

	interCmd.Flags().StringP("out-prefix", "o", "-", `out file prefix ("-" for stdout)`)
	interCmd.Flags().IntP("esti-kmer-num", "n", 100000000, "estimated kmer num length (for initializing Bloom Filter)")

}
