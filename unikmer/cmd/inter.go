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
	Short: "intersection of multiple binary files",
	Long: `intersection of multiple binary files

Tips:
  1. for comparing TWO files with really huge number of k-mers,
     you can use 'unikmer sort -m 1G -u'.

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
		sortKmers := getFlagBool(cmd, "sort")

		m := make(map[uint64]bool, mapInitSize)

		var infh *bufio.Reader
		var r *os.File
		var reader *unikmer.Reader
		var k int = -1
		var canonical bool
		var firstFile = true
		var hasInter = true
		var code uint64
		var ok bool
		var flag int
		var nfiles = len(files)
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

				if len(files) == 1 {
					reader, err = unikmer.NewReader(infh)
					checkError(err)

					k = reader.K

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
					var m2 []uint64

					if sortKmers {
						m2 = make([]uint64, 0, mapInitSize)
					} else {
						var mode uint32
						if sortKmers {
							mode |= unikmer.UNIK_SORTED
						} else if opt.Compact {
							mode |= unikmer.UNIK_COMPACT
						}
						if reader.IsCanonical() {
							mode |= unikmer.UNIK_CANONICAL
						}
						writer, err = unikmer.NewWriter(outfh, reader.K, mode)
						checkError(err)
					}

					m := make(map[uint64]struct{}, mapInitSize)
					for {
						code, err = reader.ReadCode()
						if err != nil {
							if err == io.EOF {
								break
							}
							checkError(err)
						}

						if _, ok = m[code]; !ok {
							m[code] = struct{}{}
							if sortKmers {
								m2 = append(m2, code)
							} else {
								writer.WriteCode(code) // not need to check err
							}
						}
					}

					if sortKmers {
						var mode uint32
						if reader.IsCanonical() {
							mode |= unikmer.UNIK_CANONICAL
						}
						mode |= unikmer.UNIK_SORTED
						writer, err = unikmer.NewWriter(outfh, reader.K, mode)
						checkError(err)

						writer.Number = int64(len(m2))

						if opt.Verbose {
							log.Infof("sorting %d k-mers", len(m2))
						}
						sort.Sort(unikmer.CodeSlice(m2))
						if opt.Verbose {
							log.Infof("done sorting")
						}

						for _, code = range m2 {
							writer.WriteCode(code)
						}
					}

					checkError(writer.Flush())
					if opt.Verbose {
						log.Infof("%d k-mers saved to %s", len(m), outFile)
					}
					return flagReturn
				}

				reader, err = unikmer.NewReader(infh)
				checkError(err)

				if k == -1 {
					k = reader.K
					canonical = reader.IsCanonical()
				} else if k != reader.K {
					checkError(fmt.Errorf("K (%d) of binary file '%s' not equal to previous K (%d)", reader.K, file, k))
				} else if (reader.IsCanonical()) != canonical {
					checkError(fmt.Errorf(`'canonical' flags not consistent, please check with "unikmer stats"`))
				}

				for {
					code, err = reader.ReadCode()
					if err != nil {
						if err == io.EOF {
							break
						}
						checkError(err)
					}

					if firstFile {
						m[code] = false
						continue
					}

					// mark seen kmer
					if _, ok = m[code]; ok {
						m[code] = true
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

		writer, err := unikmer.NewWriter(outfh, k, mode)
		checkError(err)

		if sortKmers {
			writer.Number = int64(len(m))
		}

		if len(m) > 0 {
			if sortKmers {
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
			} else {
				for code = range m {
					// not need to check err
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
