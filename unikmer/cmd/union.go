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

package cmd

import (
	"bufio"
	"fmt"
	"io"
	"os"
	"strings"

	"github.com/pkg/errors"
	"github.com/shenwei356/bio/taxdump"
	"github.com/shenwei356/unik/v5"

	"github.com/spf13/cobra"
	"github.com/twotwotwo/sorts/sortutil"
)

var unionCmd = &cobra.Command{
	Use:   "union",
	Short: "Union of k-mers in multiple binary files",
	Long: `Union of k-mers in multiple binary files

Attentions:
  1. The 'canonical/scaled/hashed' flags of all files should be consistent.
  2. Input files should ALL have or don't have taxid information.

Tips:
  1. 'unikmer sort -u' is slightly faster in cost of more memory usage.
  2. For really huge number of k-mers, you can use 'unikmer sort -m 100M -u'.
  3. For large number of sorted .unik files, you can use 'unikmer merge'.

`,
	Run: func(cmd *cobra.Command, args []string) {
		opt := getOptions(cmd)

		var err error

		if opt.Verbose {
			log.Info("checking input files ...")
		}
		files := getFileListFromArgsAndFile(cmd, args, true, "infile-list", !opt.SkipFileCheck)
		if opt.Verbose {
			if len(files) == 1 && isStdin(files[0]) {
				log.Info("no files given, reading from stdin")
			} else {
				log.Infof("%d input file(s) given", len(files))
			}
		}

		checkFileSuffix(opt, extDataFile, files...)

		outFile := getFlagString(cmd, "out-prefix")
		sortKmers := getFlagBool(cmd, "sort")

		var m map[uint64]struct{}
		var taxondb *taxdump.Taxonomy
		var mt map[uint64]uint32

		if !isStdout(outFile) && !strings.HasSuffix(outFile, extDataFile) {
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

		var writer *unik.Writer

		var infh *bufio.Reader
		var r *os.File

		if len(files) == 1 {
			if opt.Verbose {
				log.Infof("directly copy the only one input file to output file")
			}
			infh, r, _, err = inStream(files[0])
			checkError(err)
			defer r.Close()

			if !isStdout(outFile) && !strings.HasSuffix(outFile, extDataFile) {
				outFile += extDataFile
			}

			_, err = io.Copy(outfh, infh)
			checkError(err)
			return
		}

		var reader0 *unik.Reader
		var code uint64
		var taxid uint32
		var lca uint32
		var k int = -1
		var canonical bool
		var hashed bool
		var hasTaxid bool
		var ok bool
		var n int
		var flag int
		var nfiles = len(files)
		for i, file := range files {
			if opt.Verbose {
				log.Infof("processing file (%d/%d): %s", i+1, nfiles, file)
			}

			flag = func() int {
				infh, r, _, err = inStream(file)
				checkError(err)
				defer r.Close()

				reader, err := unik.NewReader(infh)
				checkError(errors.Wrap(err, file))

				if k == -1 {
					reader0 = reader
					k = reader.K
					canonical = reader.IsCanonical()
					hashed = reader.IsHashed()
					hasTaxid = !opt.IgnoreTaxid && reader.HasTaxidInfo()
					if hasTaxid {
						if opt.Verbose {
							log.Infof("taxids found in file: %s", file)
						}
						mt = make(map[uint64]uint32, mapInitSize)
						taxondb = loadTaxonomy(opt, false)
					} else {
						m = make(map[uint64]struct{}, mapInitSize)
					}

					if !hasTaxid && !sortKmers {
						var mode uint32
						if sortKmers {
							mode |= unik.UnikSorted
						} else if opt.Compact && !hashed {
							mode |= unik.UnikCompact
						}
						if canonical {
							mode |= unik.UnikCanonical
						}
						if hasTaxid {
							mode |= unik.UnikIncludeTaxID
						}
						if hashed {
							mode |= unik.UnikHashed
						}
						writer, err = unik.NewWriter(outfh, k, mode)
						checkError(errors.Wrap(err, outFile))
						writer.SetMaxTaxid(opt.MaxTaxid)
					}
				} else {
					checkCompatibility(reader0, reader, file)
					if !opt.IgnoreTaxid && reader.HasTaxidInfo() != hasTaxid {
						if reader.HasTaxidInfo() {
							checkError(fmt.Errorf(`taxid information not found in previous files, but found in this: %s`, file))
						} else {
							checkError(fmt.Errorf(`taxid information found in previous files, but missing in this: %s`, file))
						}
					}
				}

				for {
					code, taxid, err = reader.ReadCodeWithTaxid()
					if err != nil {
						if err == io.EOF {
							break
						}
						checkError(errors.Wrap(err, file))
					}

					if hasTaxid {
						if lca, ok = mt[code]; !ok {
							mt[code] = taxid
						} else {
							mt[code] = taxondb.LCA(lca, taxid) // update with LCA
						}
						continue
					}

					if _, ok = m[code]; !ok {
						m[code] = struct{}{}
						n++
					}
				}

				return flagContinue
			}()

			if flag == flagReturn {
				return
			} else if flag == flagBreak {
				break
			}
		}

		if sortKmers || hasTaxid {
			var mode uint32
			if sortKmers {
				mode |= unik.UnikSorted
			} else if opt.Compact && !hashed {
				mode |= unik.UnikCompact
			}
			if canonical {
				mode |= unik.UnikCanonical
			}
			if hasTaxid {
				mode |= unik.UnikIncludeTaxID
			}
			if hashed {
				mode |= unik.UnikHashed
			}
			writer, err = unik.NewWriter(outfh, k, mode)
			checkError(err)
			writer.SetMaxTaxid(opt.MaxTaxid)

			if hasTaxid {
				n = len(mt)
			} else {
				n = len(m)
			}
			writer.Number = uint64(n)
		}

		if !sortKmers {
			if hasTaxid {
				for code, taxid = range mt {
					writer.WriteCodeWithTaxid(code, taxid)
				}
				n = len(mt)
			} else {
				for code = range m {
					writer.WriteCode(code)
				}
				n = len(m)
			}
		} else {
			if hasTaxid {
				codes := make([]uint64, len(mt))

				i := 0
				for code = range mt {
					codes[i] = code
					i++
				}

				if opt.Verbose {
					log.Infof("sorting %d k-mers", len(codes))
				}
				// sort.Sort(unikmer.CodeSlice(codes))
				sortutil.Uint64s(codes)
				if opt.Verbose {
					log.Infof("done sorting")
				}
				for _, code = range codes {
					writer.WriteCodeWithTaxid(code, mt[code])
				}
				n = len(mt)
			} else {
				codes := make([]uint64, len(m))

				i := 0
				for code = range m {
					codes[i] = code
					i++
				}

				if opt.Verbose {
					log.Infof("sorting %d k-mers", len(codes))
				}
				// sort.Sort(unikmer.CodeSlice(codes))
				sortutil.Uint64s(codes)
				if opt.Verbose {
					log.Infof("done sorting")
				}

				for _, code = range codes {
					writer.WriteCode(code)
				}
				n = len(m)
			}
		}

		checkError(writer.Flush())
		if opt.Verbose {
			log.Infof("%d k-mers saved to %s", n, outFile)
		}
	},
}

func init() {
	RootCmd.AddCommand(unionCmd)

	unionCmd.Flags().StringP("out-prefix", "o", "-", `out file prefix ("-" for stdout)`)
	unionCmd.Flags().BoolP("sort", "s", false, helpSort)
}
