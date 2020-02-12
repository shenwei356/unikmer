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
	"sort"

	"github.com/shenwei356/unikmer"
	"github.com/spf13/cobra"
)

// interCmd represents
var interCmd = &cobra.Command{
	Use:   "inter",
	Short: "Intersection of multiple binary files",
	Long: `Intersection of multiple binary files

Attentions:
  1. The 'canonical' flags of all files should be consistent.
  2. Input files should ALL have or don't have taxid information.
  
Tips:
  1. For comparing TWO files with really huge number of k-mers,
     you can use 'unikmer sort -u -m 100M' for each file,
     and then 'unikmer merge -' from them.

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
		var nfiles = len(files)

		if nfiles <= 1 {
			checkError(fmt.Errorf("two or more files needed"))
		}

		outFile := getFlagString(cmd, "out-prefix")
		sortKmers := getFlagBool(cmd, "sort")

		var m map[uint64]bool
		var taxondb *unikmer.Taxonomy
		var mt map[uint64]uint32

		m = make(map[uint64]bool, mapInitSize) // if has taxid, this map is used for marking common elements

		var infh *bufio.Reader
		var r *os.File
		var reader *unikmer.Reader
		var k int = -1
		var canonical bool
		var hasTaxid bool
		var firstFile = true
		var hasInter = true
		var code uint64
		var taxid uint32
		var ok bool
		var flag int

		for i, file := range files {
			if !firstFile && file == files[0] {
				continue
			}

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

					if hasTaxid {
						if opt.Verbose {
							log.Infof("taxids found in file: %s", file)
						}
						mt = make(map[uint64]uint32, mapInitSize)
						taxondb = loadTaxonomy(opt)
					}
				} else {
					if k != reader.K {
						checkError(fmt.Errorf("K (%d) of binary file '%s' not equal to previous K (%d)", reader.K, file, k))
					}
					if reader.IsCanonical() != canonical {
						checkError(fmt.Errorf(`'canonical' flags not consistent, please check with "unikmer stats"`))
					}
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
						checkError(err)
					}

					if firstFile {
						m[code] = false
						if hasTaxid {
							mt[code] = taxid
						}
						continue
					}

					// mark seen kmer
					if _, ok = m[code]; ok {
						m[code] = true
						if hasTaxid {
							mt[code] = taxondb.LCA(mt[code], taxid)
						}
					}
				}

				if firstFile {
					firstFile = false
					return flagContinue
				}

				// remove unseen kmers
				for code = range m {
					if m[code] {
						m[code] = false
					} else {
						delete(m, code)
						if hasTaxid {
							delete(mt, code)
						}
					}
				}

				if opt.Verbose {
					log.Infof("%d k-mers remain", len(m))
				}
				if len(m) == 0 {
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
		if sortKmers {
			mode |= unikmer.UNIK_SORTED
		} else if opt.Compact {
			mode |= unikmer.UNIK_COMPACT
		}
		if canonical {
			mode |= unikmer.UNIK_CANONICAL
		}
		if hasTaxid {
			mode |= unikmer.UNIK_INCLUDETAXID
		}

		writer, err := unikmer.NewWriter(outfh, k, mode)
		checkError(err)
		writer.SetMaxTaxid(opt.MaxTaxid) // follow taxondb

		if sortKmers {
			writer.Number = int64(len(m))
		}

		if !sortKmers {
			if hasTaxid {
				for code, taxid = range mt {
					writer.WriteCodeWithTaxid(code, taxid)
				}
			} else {
				for code = range m {
					writer.WriteCode(code)
				}
			}
		} else if len(m) > 0 {
			if hasTaxid {
				codesTaxids := make([]unikmer.CodeTaxid, len(mt))

				i := 0
				for code, taxid := range mt {
					codesTaxids[i] = unikmer.CodeTaxid{Code: code, Taxid: taxid}
					i++
				}

				if opt.Verbose {
					log.Infof("sorting %d k-mers", len(codesTaxids))
				}
				sort.Sort(unikmer.CodeTaxidSlice(codesTaxids))
				if opt.Verbose {
					log.Infof("done sorting")
				}

				for _, codeT := range codesTaxids {
					writer.WriteCodeWithTaxid(codeT.Code, codeT.Taxid)
				}
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
				sort.Sort(unikmer.CodeSlice(codes))
				if opt.Verbose {
					log.Infof("done sorting")
				}

				for _, code = range codes {
					writer.WriteCode(code)
				}
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
	interCmd.Flags().BoolP("sort", "s", false, helpSort)
}
