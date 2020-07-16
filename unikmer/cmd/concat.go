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

// concatCmd represents
var concatCmd = &cobra.Command{
	Use:   "concat",
	Short: "Concatenate multiple binary files without removing duplicates",
	Long: `Concatenate multiple binary files without removing duplicates

Attentions:
  1. The 'canonical' flags of all files should be consistent.
  2. Input files should ALL have or don't have taxid information.

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

		checkFileSuffix(extDataFile, files...)

		outFile := getFlagString(cmd, "out-prefix")
		sortedKmers := getFlagBool(cmd, "sorted")
		globalTaxid := getFlagUint32(cmd, "taxid")
		hasGlobalTaxid := globalTaxid > 0

		if hasGlobalTaxid && opt.Verbose {
			log.Warningf("discarding all taxids and assigning new global taxid: %d", globalTaxid)
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

		var writer *unikmer.Writer

		var infh *bufio.Reader
		var r *os.File
		var reader *unikmer.Reader
		var code uint64
		var taxid uint32
		var k int = -1
		var canonical bool
		var hasTaxid bool
		var flag int
		var n int64
		var nfiles = len(files)
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
					canonical = reader.IsCanonical()
					hasTaxid = !opt.IgnoreTaxid && reader.HasTaxidInfo()

					var mode uint32
					if sortedKmers {
						mode |= unikmer.UNIK_SORTED
					} else if opt.Compact {
						mode |= unikmer.UNIK_COMPACT
					}
					if canonical {
						mode |= unikmer.UNIK_CANONICAL
					}
					if hasTaxid && !hasGlobalTaxid {
						mode |= unikmer.UNIK_INCLUDETAXID
					}
					writer, err = unikmer.NewWriter(outfh, k, mode)
					checkError(err)
					writer.SetMaxTaxid(maxUint32N(reader.GetTaxidBytesLength())) // follow reader
					if hasGlobalTaxid {
						writer.SetGlobalTaxid(globalTaxid)
					}
				} else {
					if k != reader.K {
						checkError(fmt.Errorf("K (%d) of binary file '%s' not equal to previous K (%d)", reader.K, file, k))
					}
					if reader.IsCanonical() != canonical {
						checkError(fmt.Errorf(`'canonical' flags not consistent, please check with "unikmer stats"`))
					}
					if !hasGlobalTaxid && !opt.IgnoreTaxid && reader.HasTaxidInfo() != hasTaxid {
						if reader.HasTaxidInfo() {
							checkError(fmt.Errorf(`taxid information not found in previous files, but found in this: %s`, file))
						} else {
							checkError(fmt.Errorf(`taxid information found in previous files, but missing in this: %s`, file))
						}

					}
				}

				if hasGlobalTaxid {
					for {
						code, _, err = reader.ReadCodeWithTaxid()
						if err != nil {
							if err == io.EOF {
								break
							}
							checkError(err)
						}

						checkError(writer.WriteCode(code))
						n++
					}

					return flagContinue
				}

				for {
					code, taxid, err = reader.ReadCodeWithTaxid()
					if err != nil {
						if err == io.EOF {
							break
						}
						checkError(err)
					}

					checkError(writer.WriteCodeWithTaxid(code, taxid))
					n++
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
	RootCmd.AddCommand(concatCmd)

	concatCmd.Flags().StringP("out-prefix", "o", "-", `out file prefix ("-" for stdout)`)
	concatCmd.Flags().BoolP("sorted", "s", false, "input k-mers are sorted")
	concatCmd.Flags().Uint32P("taxid", "t", 0, "global taxid")
}
