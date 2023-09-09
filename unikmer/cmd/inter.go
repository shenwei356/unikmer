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

	"github.com/pkg/errors"
	"github.com/shenwei356/bio/taxdump"
	"github.com/shenwei356/unik/v5"

	"github.com/spf13/cobra"
)

var interCmd = &cobra.Command{
	Use:   "inter",
	Short: "Intersection of k-mers in multiple binary files",
	Long: `Intersection of k-mers in multiple binary files

Attentions:
  0. All input files should be sorted, and output file is sorted.
  1. The 'canonical/scaled/hashed' flags of all files should be consistent.
  2. Taxid information could be inconsistent when using flag --mix-taxid.
  
Tips:
  1. For comparing TWO files with really huge number of k-mers,
     you can use 'unikmer sort -u -m 100M' for each file,
     and then 'unikmer merge -' from them.
  2. Put the smallest file in the beginning to reduce memory usage.

`,
	Run: func(cmd *cobra.Command, args []string) {
		opt := getOptions(cmd)

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

		checkFileSuffix(opt, extDataFile, files...)
		var nfiles = len(files)

		outFile := getFlagString(cmd, "out-prefix")
		mixTaxid := getFlagBool(cmd, "mix-taxid")
		var hasMixTaxid bool

		var taxondb *taxdump.Taxonomy

		mc := make([]CodeTaxid, 0, mapInitSize)
		m := make([]bool, 0, mapInitSize) // marking common elements

		var infh *bufio.Reader
		var r *os.File
		var reader0 *unik.Reader
		var k int = -1
		var canonical bool
		var hashed bool
		var hasTaxid bool
		var firstFile = true
		var hasInter = true
		var code uint64
		var taxid uint32
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

				reader, err := unik.NewReader(infh)
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

		var reader *unik.Reader
		for i, file := range files {
			if opt.Verbose {
				log.Infof("processing file (%d/%d): %s", i+1, nfiles, file)
			}

			flag = func() int {
				infh, r, _, err = inStream(file)
				checkError(err)
				defer r.Close()

				reader, err = unik.NewReader(infh)
				checkError(errors.Wrap(err, file))

				if firstFile {
					for {
						code, taxid, err = reader.ReadCodeWithTaxid()
						if err != nil {
							if err == io.EOF {
								break
							}
							checkError(errors.Wrap(err, file))
						}

						mc = append(mc, CodeTaxid{Code: code, Taxid: taxid})
						m = append(m, false)
					}
					firstFile = false
					return flagContinue
				}

				var qCode, code uint64
				var qtaxid, taxid uint32
				ii := 0
				qCode = mc[ii].Code
				qtaxid = mc[ii].Taxid

				code, taxid, err = reader.ReadCodeWithTaxid()
				if err != nil {
					if err == io.EOF {
						return flagBreak
					}
					checkError(errors.Wrap(err, file))
				}

				n := 0
				for {
					if qCode < code {
						ii++
						if ii >= len(mc) {
							break
						}
						qCode = mc[ii].Code
						qtaxid = mc[ii].Taxid
					} else if qCode == code {
						if hasMixTaxid {
							if qtaxid == 0 {
								mc[ii].Taxid = taxid
							} else if taxid == 0 {
								mc[ii].Taxid = qtaxid
							} else {
								mc[ii].Taxid = taxondb.LCA(qtaxid, taxid)
							}
						} else if hasTaxid {
							mc[ii].Taxid = taxondb.LCA(qtaxid, taxid)
						}

						m[ii] = true
						n++

						ii++
						if ii >= len(mc) {
							break
						}
						qCode = mc[ii].Code
						qtaxid = mc[ii].Taxid

						code, taxid, err = reader.ReadCodeWithTaxid()
						if err != nil {
							if err == io.EOF {
								break
							}
							checkError(errors.Wrap(err, file))
						}
					} else {
						code, taxid, err = reader.ReadCodeWithTaxid()
						if err != nil {
							if err == io.EOF {
								break
							}
							checkError(errors.Wrap(err, file))
						}
					}
				}

				mc1 := make([]CodeTaxid, 0, n)
				n = 0
				for ii, found := range m {
					if found {
						mc1 = append(mc1, mc[ii])
						n++
					}
				}
				mc = mc1
				m = make([]bool, n)

				if opt.Verbose {
					log.Infof("%d k-mers remain", n)
				}
				if n == 0 {
					hasInter = false
					return flagBreak
				}

				return flagContinue
			}()

			if flag == flagReturn {
				return
			} else if flag == flagBreak {
				break
			}
		}

		if !hasInter {
			if opt.Verbose {
				log.Infof("no intersection found")
			}
			// return
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
		mode |= unik.UnikSorted
		if canonical {
			mode |= unik.UnikCanonical
		}
		if hasTaxid || hasMixTaxid {
			mode |= unik.UnikIncludeTaxID
		}
		if hashed {
			mode |= unik.UnikHashed
		}

		writer, err := unik.NewWriter(outfh, k, mode)
		checkError(errors.Wrap(err, outFile))
		writer.SetMaxTaxid(opt.MaxTaxid) // follow taxondb

		writer.Number = uint64(len(mc))

		if hasTaxid || hasMixTaxid {
			for _, ct := range mc {
				writer.WriteCodeWithTaxid(ct.Code, ct.Taxid)
			}
		} else {
			for _, ct := range mc {
				writer.WriteCode(ct.Code)
			}
		}

		checkError(writer.Flush())
		if opt.Verbose {
			log.Infof("%d k-mers saved to %s", len(m), outFile)
		}
	},
}

func init() {
	RootCmd.AddCommand(interCmd)

	interCmd.Flags().StringP("out-prefix", "o", "-", `out file prefix ("-" for stdout)`)
	interCmd.Flags().BoolP("mix-taxid", "m", false, `allow part of files being whithout taxids`)
}
