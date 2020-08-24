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
	"path/filepath"
	"runtime"
	"sort"
	"sync"

	"github.com/pkg/errors"
	"github.com/shenwei356/util/pathutil"

	"github.com/shenwei356/unikmer"
	"github.com/spf13/cobra"
)

var splitCmd = &cobra.Command{
	Use:   "split",
	Short: "Split k-mers into sorted chunk files",
	Long: `Split k-mers into sorted chunk files

Attentions:
  1. The 'canonical' flags of all files should be consistent.
  2. Input files should ALL have or don't have taxid information.
  
Tips:
  1. You can use '-m/--chunk-size' to limit memory usage, and chunk file size
     depends on k-mers and file save mode (sorted/compact/normal).
  2. Increasing value of -j/--threads can accelerates splitting stage,
     in cost of more memory occupation.
  3. For sorted input files, the memory usage is very low and speed is fast.

`,
	Run: func(cmd *cobra.Command, args []string) {
		opt := getOptions(cmd)
		runtime.GOMAXPROCS(opt.NumCPUs)

		outDir := getFlagString(cmd, "out-dir")
		force := getFlagBool(cmd, "force")
		unique := getFlagBool(cmd, "unique")
		repeated := getFlagBool(cmd, "repeated")

		maxElem, err := ParseByteSize(getFlagString(cmd, "chunk-size"))
		if err != nil {
			checkError(fmt.Errorf("parsing byte size: %s", err))
		}
		limitMem := maxElem > 0

		var listInitSize int
		if limitMem {
			listInitSize = maxElem
		} else {
			listInitSize = mapInitSize
		}

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

		var m []uint64
		var taxondb *unikmer.Taxonomy
		var mt []unikmer.CodeTaxid

		if outDir == "" {
			if isStdin(files[0]) {
				outDir = "stdin.split"
			} else {
				outDir = files[0] + ".split"
			}
		}
		pwd, _ := os.Getwd()
		if outDir != "./" && outDir != "." && pwd != filepath.Clean(outDir) {
			existed, err := pathutil.DirExists(outDir)
			checkError(errors.Wrap(err, outDir))
			if existed {
				empty, err := pathutil.IsEmpty(outDir)
				checkError(errors.Wrap(err, outDir))
				if !empty {
					if force {
						checkError(os.RemoveAll(outDir))
					} else {
						checkError(fmt.Errorf("outDir not empty: %s, use --force to overwrite", outDir))
					}
				} else {
					checkError(os.RemoveAll(outDir))
				}
			}
			checkError(os.MkdirAll(outDir, 0777))
		}

		var infh *bufio.Reader
		var r *os.File
		var reader0 *unikmer.Reader
		var code uint64
		var taxid uint32
		var k int = -1
		var canonical bool
		var hashed bool
		var hasTaxid bool
		var mode uint32
		var flag int
		var nfiles = len(files)
		var doNotNeedSorting = false // only for ONE sorted input file

		var iTmpFile int

		var wg sync.WaitGroup
		tokens := make(chan int, opt.NumCPUs)

		// just for counting total k-mers
		chN := make(chan int64, opt.NumCPUs)
		var N int64 = 0
		done := make(chan int)
		go func() {
			for n := range chN {
				N += n

				runtime.GC()
			}
			done <- 1
		}()

		var outFile2 string
		var n int
		var writer *unikmer.Writer
		var outfh *bufio.Writer
		var gw io.WriteCloser
		var w *os.File

		for i, file := range files {

			if opt.Verbose {
				log.Infof("processing file (%d/%d): %s", i+1, nfiles, file)
			}

			flag = func() int {
				infh, r, _, err = inStream(file)
				checkError(err)
				defer r.Close()

				reader, err := unikmer.NewReader(infh)
				checkError(errors.Wrap(err, file))

				if k == -1 {
					reader0 = reader
					k = reader.K
					canonical = reader.IsCanonical()
					hashed = reader.IsHashed()
					hasTaxid = !opt.IgnoreTaxid && reader.HasTaxidInfo()
					if nfiles == 1 && reader.IsSorted() {
						doNotNeedSorting = true
					}

					if !doNotNeedSorting {
						if hasTaxid {
							if opt.Verbose {
								log.Infof("taxids found in file: %s", file)
							}
							mt = make([]unikmer.CodeTaxid, 0, listInitSize)
							taxondb = loadTaxonomy(opt, false)
						} else {
							m = make([]uint64, 0, listInitSize)
						}
					} else {
						log.Infof("sorting is not needed for ONE input file")
					}

					if hasTaxid {
						mode |= unikmer.UNIK_INCLUDETAXID
					}
					if hashed {
						mode |= unikmer.UNIK_HASHED
					}
					mode |= unikmer.UNIK_SORTED

					if doNotNeedSorting {
						iTmpFile++
						outFile2 = chunkFileName(outDir, iTmpFile)
						outfh, gw, w, err = outStream(outFile2, opt.Compress, opt.CompressionLevel)
						checkError(err)

						writer, err = unikmer.NewWriter(outfh, k, mode)
						checkError(errors.Wrap(err, outFile2))
						writer.SetMaxTaxid(maxUint32N(reader.GetTaxidBytesLength())) // follow reader
						if opt.Verbose {
							log.Infof("[chunk %d] begin writing k-mers to: %s", iTmpFile, outFile2)
						}
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

					if doNotNeedSorting {
						writer.WriteCodeWithTaxid(code, taxid)
						n++

						if limitMem && n >= maxElem {
							if opt.Verbose {
								log.Infof("[chunk %d] %d k-mers saved to: %s", iTmpFile, n, outFile2)
							}
							chN <- int64(n)

							writer.Flush()
							outfh.Flush()
							if gw != nil {
								gw.Close()
							}
							w.Close()

							iTmpFile++
							outFile2 = chunkFileName(outDir, iTmpFile)
							outfh, gw, w, err = outStream(outFile2, opt.Compress, opt.CompressionLevel)
							checkError(err)

							writer, err = unikmer.NewWriter(outfh, k, mode)
							checkError(errors.Wrap(err, outFile2))
							writer.SetMaxTaxid(maxUint32N(reader.GetTaxidBytesLength())) // follow reader

							if opt.Verbose {
								log.Infof("[chunk %d] begin writing k-mers to: %s", iTmpFile, outFile2)
							}

							n = 0
						}

						continue
					}

					if hasTaxid {
						mt = append(mt, unikmer.CodeTaxid{Code: code, Taxid: taxid})
					} else {
						m = append(m, code)
					}

					if limitMem && (len(m) >= maxElem || len(mt) >= maxElem) {
						iTmpFile++
						outFile1 := chunkFileName(outDir, iTmpFile)

						wg.Add(1)
						tokens <- 1
						go func(m []uint64, mt []unikmer.CodeTaxid, iTmpFile int, outFile string) {
							defer func() {
								wg.Done()
								<-tokens
							}()

							if hasTaxid {
								if opt.Verbose {
									log.Infof("[chunk %d] sorting %d k-mers", iTmpFile, len(mt))
								}
								sort.Sort(unikmer.CodeTaxidSlice(mt))
							} else {
								if opt.Verbose {
									log.Infof("[chunk %d] sorting %d k-mers", iTmpFile, len(m))
								}
								sort.Sort(unikmer.CodeSlice(m))
							}

							var _n int64
							if hasTaxid {
								_n = dumpCodesTaxids2File(mt, taxondb, k, mode, outFile, opt, unique, repeated)
							} else {
								_n = dumpCodes2File(m, k, mode, outFile, opt, unique, repeated)
							}
							if opt.Verbose {
								log.Infof("[chunk %d] %d k-mers saved to %s", iTmpFile, _n, outFile)
							}
							chN <- int64(len(m))
						}(m, mt, iTmpFile, outFile1)

						if hasTaxid {
							mt = make([]unikmer.CodeTaxid, 0, listInitSize)
						} else {
							m = make([]uint64, 0, listInitSize)
						}

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

		// dump remaining k-mers to file

		if doNotNeedSorting {
			writer.Flush()
			outfh.Flush()
			if gw != nil {
				gw.Close()
			}
			w.Close()

			if n > 0 {
				if opt.Verbose {
					log.Infof("[chunk %d] %d k-mers saved to %s", iTmpFile, n, outFile2)
				}
				chN <- int64(n)
			} else {
				os.Remove(outFile2)
				iTmpFile--
			}
		}

		if len(m) > 0 || len(mt) > 0 {
			iTmpFile++
			outFile1 := chunkFileName(outDir, iTmpFile)

			wg.Add(1)
			tokens <- 1
			go func(m []uint64, mt []unikmer.CodeTaxid, iTmpFile int, outFile string) {
				defer func() {
					wg.Done()
					<-tokens
				}()

				if opt.Verbose {
					log.Infof("[chunk %d] sorting %d k-mers", iTmpFile, len(m))
				}
				sort.Sort(unikmer.CodeSlice(m))
				if opt.Verbose {
					log.Infof("[chunk %d] done sorting", iTmpFile)
					log.Infof("[chunk %d] writing to file: %s", iTmpFile, outFile)
				}

				var _n int64
				if hasTaxid {
					_n = dumpCodesTaxids2File(mt, taxondb, k, mode, outFile, opt, unique, repeated)
				} else {
					_n = dumpCodes2File(m, k, mode, outFile, opt, unique, repeated)
				}
				if opt.Verbose {
					log.Infof("[chunk %d] %d k-mers saved to %s", iTmpFile, _n, outFile)
				}
				chN <- int64(len(m))

			}(m, mt, iTmpFile, outFile1)
		}

		// wait all k-mers being wrote to files
		wg.Wait()
		close(chN)
		<-done

		if opt.Verbose {
			log.Infof("%d chunk files with total %d k-mers saved to dir: %s", iTmpFile, N, outDir)
		}
	},
}

func init() {
	RootCmd.AddCommand(splitCmd)

	splitCmd.Flags().StringP("out-dir", "O", "", `output directory`)
	splitCmd.Flags().StringP("chunk-size", "m", "", `split input into chunks of N k-mers, supports K/M/G suffix, type "unikmer sort -h" for detail`)
	splitCmd.Flags().BoolP("force", "", false, `overwrite output directory`)
	splitCmd.Flags().BoolP("unique", "u", false, `split for further removing duplicated k-mers`)
	splitCmd.Flags().BoolP("repeated", "d", false, `split for further printing duplicate k-mers`)
}
