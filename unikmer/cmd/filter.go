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

// filterCmd represents
var filterCmd = &cobra.Command{
	Use:   "filter",
	Short: "Filter low-complexity k-mers",
	Long: `Filter low-complexity k-mers (experimental)

Attentions:
  1. This command only detects single base repeat now.

`,
	Run: func(cmd *cobra.Command, args []string) {
		opt := getOptions(cmd)
		runtime.GOMAXPROCS(opt.NumCPUs)

		var err error

		if opt.Verbose {
			log.Info("checking input files ...")
		}
		files := getFileListFromArgsAndFile(cmd, args, true, "infile-list", true)
		if opt.Verbose {
			if len(files) == 1 && isStdin(files[0]) {
				log.Info("no files given, reading from stdin")
			} else {
				log.Infof("%d input file(s) given", len(files))
			}
		}

		if len(files) > 1 {
			checkError(fmt.Errorf("no more than one file should be given"))
		}

		checkFileSuffix(extDataFile, files...)

		outFile := getFlagString(cmd, "out-prefix")
		threshold := getFlagNonNegativeInt(cmd, "threshold")
		invert := getFlagBool(cmd, "invert")
		window := getFlagPositiveInt(cmd, "window")
		penaltyS := getFlagInt(cmd, "penalty-s")
		penaltyD := getFlagInt(cmd, "penalty-d")

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

		var writer *unikmer.Writer

		var infh *bufio.Reader
		var r *os.File
		var reader *unikmer.Reader
		var code uint64
		var taxid uint32
		var k int = -1
		var canonical bool
		var flag int
		var nfiles = len(files)
		var hit bool
		var n int64
		var scores []int
		for i, file := range files {
			if opt.Verbose {
				log.Infof("processing file (%d/%d): %s", i+1, nfiles, file)
			}

			flag = func() int {
				infh, r, _, err = inStream(file)
				checkError(err)
				defer r.Close()

				reader, err = unikmer.NewReader(infh)
				checkError(err)

				if k == -1 {
					k = reader.K
					if window > k {
						log.Warningf("window size (%d) is bigger than k (%d)", window, k)
						window = k
					}
					canonical = reader.IsCanonical()

					scores = make([]int, k)

					writer, err = unikmer.NewWriter(outfh, k, reader.Flag)
					checkError(err)
					writer.SetMaxTaxid(maxUint32N(reader.GetTaxidBytesLength())) // follow reader
				} else {
					if k != reader.K {
						checkError(fmt.Errorf("K (%d) of binary file '%s' not equal to previous K (%d)", reader.K, file, k))
					}
					if reader.IsCanonical() != canonical {
						checkError(fmt.Errorf(`'canonical' flags not consistent, please check with "unikmer stats"`))
					}
				}

				for {
					code, taxid, err = reader.ReadCodeWithTaxid()
					if err != nil {
						if err == io.EOF {
							break
						}
						checkError(err)
					}

					hit = filterCode(code, k, penaltyS, penaltyD, threshold, window, &scores)

					if invert {
						if !hit {
							continue
						}
					} else if hit {
						continue
					}

					n++
					writer.WriteCodeWithTaxid(code, taxid) // not need to check err
				}

				return flagContinue
			}()

			if flag == flagReturn {
				return
			} else if flag == flagBreak {
				break
			}
		}

		checkError(writer.Flush())
		if opt.Verbose {
			log.Infof("%d k-mers saved to %s", n, outFile)
		}
	},
}

func init() {
	RootCmd.AddCommand(filterCmd)

	filterCmd.Flags().StringP("out-prefix", "o", "-", `out file prefix ("-" for stdout)`)
	filterCmd.Flags().BoolP("invert", "v", false, `invert result, i.e., output low-complexity k-mers`)
	filterCmd.Flags().IntP("threshold", "t", 15, `penalty threshold for filter, higher is stricter`)
	filterCmd.Flags().IntP("window", "w", 7, `window size for checking penalty`)
	filterCmd.Flags().IntP("penalty-s", "s", 3, `penalty for successive bases`)
	filterCmd.Flags().IntP("penalty-d", "d", 1, `penalty for different bases`)
}

func filterCode(code uint64, k int, penaltyS int, penaltyD int, threshold int, window int, scores *[]int) bool {
	// compute scores
	var last, c uint64
	last = 356
	for i := 0; i < k; i++ {
		c = code & 3
		if i > 0 {
			if c == last {
				(*scores)[i] = penaltyS // successive
			} else {
				(*scores)[i] = penaltyD // different
			}
		} else {
			(*scores)[i] = penaltyD
		}
		last = c
		code >>= 2
	}

	// check score in sliding window
	var s, pre int
	iLast := k - window - 1
	if iLast < 0 {
		iLast = 0
	}
	for i := 0; i <= iLast; i++ {
		if i == 0 {
			for j := 0; j < window; j++ {
				s += (*scores)[j]
			}
		} else { // update score
			s = s - pre + (*scores)[i+window-1]
		}
		pre = (*scores)[i]
		if s >= threshold {
			return true
		}
	}

	return false
}
