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
	"runtime"

	"github.com/pkg/errors"
	"github.com/shenwei356/unikmer"
	"github.com/spf13/cobra"
	"github.com/twotwotwo/sorts"
	"github.com/twotwotwo/sorts/sortutil"
)

var commonCmd = &cobra.Command{
	Use:   "common",
	Short: "Find k-mers shared by most of multiple binary files",
	Long: `Find k-mers shared by most of multiple binary files

This command is similar to "unikmer inter" but with looser restriction,
k-mers shared by some number/proportion of multiple files are outputted.

Attentions:
  0. All input files should be sorted, and output file is sorted.
  1. The 'canonical/scaled/hashed' flags of all files should be consistent.
  2. Taxid information could be inconsistent when using flag --mix-taxid.
  3. At most 65535 input files allowed.
  
Tips:
  1. For comparing TWO files with really huge number of k-mers,
     you can use 'unikmer sort -u -m 100M' for each file,
     and then 'unikmer merge -' from them.
  2. Put the smallest file in the beginning to reduce memory usage.

`,
	Run: func(cmd *cobra.Command, args []string) {
		opt := getOptions(cmd)
		runtime.GOMAXPROCS(opt.NumCPUs)
		sorts.MaxProcs = opt.NumCPUs

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
		if len(files) > 65535 {
			checkError(fmt.Errorf("at most 65535 input files allowed"))
		}

		checkFileSuffix(opt, extDataFile, files...)
		var nfiles = len(files)

		outFile := getFlagString(cmd, "out-prefix")
		mixTaxid := getFlagBool(cmd, "mix-taxid")
		var hasMixTaxid bool

		proportion := getFlagFloat64(cmd, "proportion")
		if proportion <= 0 || proportion > 1 {
			checkError(fmt.Errorf("value of -p/--proportion should be in range of (0, 1]"))
		}

		number := getFlagNonNegativeInt(cmd, "number")

		var threshold uint16

		if number == 0 {
			threshold = uint16(float64(len(files)) * proportion)
			if opt.Verbose {
				log.Infof("searching k-mers shared by >= %d (%f) files ...", threshold, proportion)
			}
		} else {
			threshold = uint16(number)
			if opt.Verbose {
				log.Infof("searching k-mers shared by >= %d files ...", threshold)
			}
		}

		var taxondb *unikmer.Taxonomy

		var mt map[uint64]uint32 // kmer -> taxid

		counts := make(map[uint64]uint16, mapInitSize) // kmer -> #files

		var infh *bufio.Reader
		var r *os.File
		var reader0 *unikmer.Reader
		var k int = -1
		var canonical bool
		var hashed bool
		var hasTaxid bool
		var firstFile = true
		var flag int

		if len(files) == 1 {
			if opt.Verbose {
				log.Infof("directly copy the only one input file to output file")
			}
			infh, r, _, err = inStream(files[0])
			checkError(err)
			defer r.Close()

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

			_, err = io.Copy(outfh, infh)
			checkError(err)
			return
		}

		// checking files
		for _, file := range files {
			if opt.NoCheckFile {
				break
			}

			if isStdin(file) {
				continue
			}
			func() {
				infh, r, _, err = inStream(file)
				checkError(err)
				defer r.Close()

				reader, err := unikmer.NewReader(infh)
				checkError(errors.Wrap(err, file))

				if !reader.IsSorted() {
					checkError(fmt.Errorf("input file should be sorted: %s", file))
				}

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
						taxondb = loadTaxonomy(opt, false)
						mt = make(map[uint64]uint32, mapInitSize)
					}
				} else {
					checkCompatibility(reader0, reader, file)
					if !opt.IgnoreTaxid && reader.HasTaxidInfo() != hasTaxid {
						if mixTaxid {
							hasMixTaxid = true
							if opt.Verbose {
								log.Infof("part of files have no taxid information")
							}
						} else if reader.HasTaxidInfo() {
							checkError(fmt.Errorf(`taxid information not found in previous files, but found in this: %s`, file))
						} else {
							checkError(fmt.Errorf(`taxid information found in previous files, but missing in this: %s`, file))
						}
					}
				}
			}()
		}

		var reader *unikmer.Reader
		for i, file := range files {
			if opt.Verbose {
				log.Infof("processing file (%d/%d): %s", i+1, nfiles, file)
			}

			flag = func() int {
				infh, r, _, err = inStream(file)
				checkError(err)
				defer r.Close()

				var code uint64
				var taxid, lca uint32
				var ok bool

				reader, err = unikmer.NewReader(infh)
				checkError(errors.Wrap(err, file))

				if firstFile {
					if hasTaxid {
						for {
							code, taxid, err = reader.ReadCodeWithTaxid()
							if err != nil {
								if err == io.EOF {
									break
								}
								checkError(errors.Wrap(err, file))
							}

							mt[code] = taxid
							counts[code] = 1
						}
					} else {
						for {
							code, taxid, err = reader.ReadCodeWithTaxid()
							if err != nil {
								if err == io.EOF {
									break
								}
								checkError(errors.Wrap(err, file))
							}

							counts[code] = 1
						}
					}

					firstFile = false
					return flagContinue
				}

				if hasTaxid {
					for {
						code, taxid, err = reader.ReadCodeWithTaxid()
						if err != nil {
							if err == io.EOF {
								break
							}
							checkError(errors.Wrap(err, file))
						}

						if lca, ok = mt[code]; !ok {
							mt[code] = taxid
						} else {
							mt[code] = taxondb.LCA(lca, taxid) // update with LCA
						}
						counts[code]++
					}
					return flagContinue
				}

				for {
					code, taxid, err = reader.ReadCodeWithTaxid()
					if err != nil {
						if err == io.EOF {
							break
						}
						checkError(errors.Wrap(err, file))
					}

					counts[code]++
				}

				return flagContinue
			}()

			if flag == flagReturn {
				return
			} else if flag == flagBreak {
				break
			}
		}

		// output

		if opt.Verbose {
			log.Infof("exporting k-mers")
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
		mode |= unikmer.UnikSorted
		if canonical {
			mode |= unikmer.UnikCanonical
		}
		if hasTaxid || hasMixTaxid {
			mode |= unikmer.UnikIncludeTaxID
		}
		if hashed {
			mode |= unikmer.UnikHashed
		}

		writer, err := unikmer.NewWriter(outfh, k, mode)
		checkError(errors.Wrap(err, outFile))
		writer.SetMaxTaxid(opt.MaxTaxid) // follow taxondb

		codes := make([]uint64, 0, mapInitSize)

		for code, count := range counts {
			if count >= threshold {
				codes = append(codes, code)
			}
		}

		writer.Number = uint64(len(codes))

		if opt.Verbose && len(codes) == 0 {
			log.Infof("no shared k-mers found")
		}

		// sort.Sort(unikmer.CodeSlice(codes))
		sortutil.Uint64s(codes)

		if hasTaxid || hasMixTaxid {
			for _, code := range codes {
				writer.WriteCodeWithTaxid(code, mt[code])
			}
		} else {
			for _, code := range codes {
				writer.WriteCode(code)
			}
		}

		checkError(writer.Flush())
		if opt.Verbose {
			log.Infof("%d k-mers saved to %s", len(codes), outFile)
		}
	},
}

func init() {
	RootCmd.AddCommand(commonCmd)

	commonCmd.Flags().StringP("out-prefix", "o", "-", `out file prefix ("-" for stdout)`)
	commonCmd.Flags().BoolP("mix-taxid", "m", false, `allow part of files being whithout taxids`)
	commonCmd.Flags().Float64P("proportion", "p", 1, `minimum proportion of files that share a k-mer`)
	commonCmd.Flags().IntP("number", "n", 0, `minimum number of files that share a k-mer (overides -p/--proportion)`)
}
