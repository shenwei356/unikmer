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

	"github.com/shenwei356/unikmer"
	"github.com/shenwei356/util/pathutil"
	"github.com/spf13/cobra"
)

// sortCmd represents
var sortCmd = &cobra.Command{
	Use:   "sort",
	Short: "Sort k-mers in binary files to reduce file size",
	Long: `Sort k-mers in binary files to reduce file size

Attentions:
  1. The 'canonical' flags of all files should be consistent.
  2. Input files should ALL have or don't have taxid information.

Notes:
  1. When sorting from large number of files, this command is equivalent to
     'unikmer split' + 'unikmer merge'.

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

		outFile0 := getFlagString(cmd, "out-prefix")
		unique := getFlagBool(cmd, "unique")
		repeated := getFlagBool(cmd, "repeated")
		tmpDir := getFlagString(cmd, "tmp-dir")
		maxOpenFiles := getFlagPositiveInt(cmd, "max-open-files")
		keepTmpDir := getFlagBool(cmd, "keep-tmp-dir")
		force := getFlagBool(cmd, "force")

		if unique && repeated {
			checkError(fmt.Errorf("flag -u/--unique overides -d/--repeated, don't provide both"))
		}

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

		checkFileSuffix(extDataFile, files...)

		var m []uint64
		var taxondb *unikmer.Taxonomy
		var mt []unikmer.CodeTaxid

		outFile := outFile0
		if !isStdout(outFile) {
			outFile += extDataFile
		}

		if limitMem {
			if isStdout(outFile0) {
				tmpDir = filepath.Join(tmpDir, "stdout.tmp")
			} else {
				tmpDir = filepath.Join(tmpDir, filepath.Base(outFile0)+".tmp")
			}

			existed, err := pathutil.DirExists(tmpDir)
			checkError(err)
			if existed {
				empty, err := pathutil.IsEmpty(tmpDir)
				checkError(err)
				if !empty {
					if force {
						checkError(os.RemoveAll(tmpDir))
					} else {
						checkError(fmt.Errorf("tmp dir not empty: %s, choose another one or use -f (--force) to overwrite", tmpDir))
					}
				} else {
					checkError(os.RemoveAll(tmpDir))
				}
			}
			checkError(os.MkdirAll(tmpDir, 0777))
		}

		var writer *unikmer.Writer

		var infh *bufio.Reader
		var r *os.File
		var reader *unikmer.Reader
		var code uint64
		var taxid uint32
		var k int = -1
		var canonical bool
		var hasTaxid bool
		var mode uint32
		var firstFile = true
		var flag int
		var nfiles = len(files)

		var tmpFiles []string
		var iTmpFile int
		var hasTmpFile bool

		var wg sync.WaitGroup
		tokens := make(chan int, opt.NumCPUs)

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
						mt = make([]unikmer.CodeTaxid, 0, listInitSize)
						taxondb = loadTaxonomy(opt)
					} else {
						m = make([]uint64, 0, listInitSize)
					}

					if canonical {
						mode |= unikmer.UNIK_CANONICAL
					}
					if hasTaxid {
						mode |= unikmer.UNIK_INCLUDETAXID
					}
					mode |= unikmer.UNIK_SORTED
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

					if hasTaxid {
						mt = append(mt, unikmer.CodeTaxid{Code: code, Taxid: taxid})
					} else {
						m = append(m, code)
					}

					if limitMem && (len(m) >= maxElem || len(mt) >= maxElem) {
						if !hasTmpFile {
							if opt.Verbose {
								log.Info()
								log.Infof("======= Stage 1: spliting k-mers into chunks =======")
							}

							tmpFiles = make([]string, 0, 10)
							hasTmpFile = true
						}
						iTmpFile++
						outFile1 := chunkFileName(tmpDir, iTmpFile)
						tmpFiles = append(tmpFiles, outFile1)

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
							if opt.Verbose {
								log.Infof("[chunk %d] done sorting", iTmpFile)
								log.Infof("[chunk %d] writing to file: %s", iTmpFile, outFile)
							}

							var _n int64
							if hasTaxid {
								_n = dumpCodesTaxids2File(mt, k, mode, outFile, opt, unique, repeated)
							} else {
								_n = dumpCodes2File(m, k, mode, outFile, opt, unique, repeated)
							}
							if opt.Verbose {
								log.Infof("[chunk %d] %d k-mers saved to tmp file: %s", iTmpFile, _n, outFile)
							}
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

		if hasTmpFile {
			// dump remaining k-mers to file
			if len(m) > 0 || len(mt) > 0 {
				iTmpFile++
				outFile1 := chunkFileName(tmpDir, iTmpFile)
				tmpFiles = append(tmpFiles, outFile1)

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
					if opt.Verbose {
						log.Infof("[chunk %d] done sorting", iTmpFile)
						log.Infof("[chunk %d] writing to file: %s", iTmpFile, outFile)
					}

					var _n int64
					if hasTaxid {
						_n = dumpCodesTaxids2File(mt, k, mode, outFile, opt, unique, repeated)
					} else {
						_n = dumpCodes2File(m, k, mode, outFile, opt, unique, repeated)
					}
					if opt.Verbose {
						log.Infof("[chunk %d] %d k-mers saved to tmp file: %s", iTmpFile, _n, outFile)
					}
				}(m, mt, iTmpFile, outFile1)
			}

			// wait all k-mers being wrote to files
			wg.Wait()

			// merge sort

			// this implemention was heavily inspired by https://github.com/oxtoacart/emsort/blob/master/emsort.go

			files = make([]string, len(tmpFiles))
			copy(files, tmpFiles)
			tmpFiles = make([]string, 0, 10)

			var n int64
			var _files []string
			if len(files) < maxOpenFiles {
				if opt.Verbose {
					log.Info()
					log.Infof("======= Stage 2: merging from %d chunks =======", len(files))
				}
				n, _ = mergeChunksFile(opt, taxondb, files, outFile, k, mode, unique, repeated, true)
			} else {
				if opt.Verbose {
					log.Info()
					log.Infof("======= Stage 2: merging from %d chunks (round: 1/2) =======", len(files))
				}

				_files = make([]string, 0, maxOpenFiles)
				for _, file := range files {
					_files = append(_files, file)
					if len(_files) == maxOpenFiles {
						iTmpFile++
						outFile1 := chunkFileName(tmpDir, iTmpFile)

						if opt.Verbose {
							log.Infof("[chunk %d] sorting k-mers from %d tmp files", iTmpFile, len(_files))
						}
						n, _ := mergeChunksFile(opt, taxondb, _files, outFile1, k, mode, unique, repeated, false)
						if opt.Verbose {
							log.Infof("%d k-mers saved to tmp file: %s", n, outFile1)
						}
						tmpFiles = append(tmpFiles, outFile1)
						_files = make([]string, 0, maxOpenFiles)
					}
				}
				if len(_files) > 0 {
					iTmpFile++
					outFile1 := chunkFileName(tmpDir, iTmpFile)

					if opt.Verbose {
						log.Infof("[chunk %d] sorting k-mers from %d tmp files", iTmpFile, len(_files))
					}
					n, _ := mergeChunksFile(opt, taxondb, _files, outFile1, k, mode, unique, repeated, false)
					if opt.Verbose {
						log.Infof("%d k-mers saved to tmp file: %s", n, outFile1)
					}
					tmpFiles = append(tmpFiles, outFile1)
				}
				if opt.Verbose {
					log.Info()
					log.Infof("======= Stage 3: merging from %d chunks (round: 2/2) =======", len(tmpFiles))
				}
				n, _ = mergeChunksFile(opt, taxondb, tmpFiles, outFile, k, mode, unique, repeated, true)
			}
			if opt.Verbose {
				log.Infof("%d k-mers saved to %s", n, outFile)
			}

			// cleanning

			if keepTmpDir {
				return
			}

			if opt.Verbose {
				log.Infof("removing %d intermediate files", len(tmpFiles)+len(files))
			}
			for _, file := range append(files, tmpFiles...) {
				err := os.Remove(file)
				if err != nil {
					checkError(fmt.Errorf("fail to remove intermediate file: %s", file))
				}
			}
			if opt.Verbose {
				log.Infof("removing tmp dir: %s", tmpDir)
			}
			err = os.Remove(tmpDir)
			if err != nil {
				checkError(fmt.Errorf("fail to remove temp directory, please manually delete it: %s", tmpDir))
			}

			return
		}

		// all k-mers are stored in memory

		if hasTaxid {
			if opt.Verbose {
				log.Infof("sorting %d k-mers", len(mt))
			}
			sort.Sort(unikmer.CodeTaxidSlice(mt))
		} else {
			if opt.Verbose {
				log.Infof("sorting %d k-mers", len(m))
			}
			sort.Sort(unikmer.CodeSlice(m))
		}
		if opt.Verbose {
			log.Infof("done sorting")
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
		writer, err = unikmer.NewWriter(outfh, k, mode)
		checkError(err)
		writer.SetMaxTaxid(opt.MaxTaxid) // follow taxondb

		var n int
		if hasTaxid {
			if unique {
				var last uint64 = ^uint64(0)
				var first bool = true
				var lca uint32
				for _, codeT := range mt {
					// same k-mer, compute LCA and handle it later
					if codeT.Code == last {
						lca = taxondb.LCA(codeT.Taxid, lca)
						continue
					}

					if first { // just ignore first code, faster than comparing code or slice index, I think
						first = false
					} else { // when meeting new k-mer, output previous one
						writer.WriteCodeWithTaxid(last, lca)
						n++
					}

					last = codeT.Code
					lca = codeT.Taxid
				}
				// do not forget the last one
				writer.WriteCodeWithTaxid(last, lca)
				n++
			} else if repeated {
				var last uint64 = ^uint64(0)
				var count int = 1
				var lca uint32
				for _, codeT := range mt {
					// same k-mer, compute LCA and handle it later
					if codeT.Code == last {
						lca = taxondb.LCA(codeT.Taxid, lca)
						count++
						continue
					}

					if count > 1 { // repeated
						writer.WriteCodeWithTaxid(last, lca)
						n++
						count = 1
					}
					last = codeT.Code
					lca = codeT.Taxid
				}
				if count > 1 { // last one
					writer.WriteCodeWithTaxid(last, lca)
					n++
					count = 0
				}
			} else {
				writer.Number = int64(len(mt))
				for _, codeT := range mt {
					writer.WriteCodeWithTaxid(codeT.Code, codeT.Taxid)
				}
				n = len(mt)
			}
		} else {
			if unique {
				var last uint64 = ^uint64(0)
				for _, code := range m {
					if code == last {
						continue
					}
					last = code
					writer.WriteCode(code)
					n++
				}
			} else if repeated {
				var last uint64 = ^uint64(0)
				var count int
				for _, code := range m {
					if code == last {
						if count == 1 { // write once
							writer.WriteCode(code)
							n++
						}
						count++
					} else {
						last = code
						count = 1
					}
				}
			} else {
				writer.Number = int64(len(m))
				for _, code := range m {
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
	RootCmd.AddCommand(sortCmd)

	sortCmd.Flags().StringP("out-prefix", "o", "-", `out file prefix ("-" for stdout)`)
	sortCmd.Flags().BoolP("unique", "u", false, `remove duplicated k-mers`)
	sortCmd.Flags().BoolP("repeated", "d", false, `only print duplicate k-mers`)
	sortCmd.Flags().StringP("chunk-size", "m", "", `split input into chunks of N k-mers, supports K/M/G suffix, type "unikmer sort -h" for detail`)
	sortCmd.Flags().StringP("tmp-dir", "t", "./", `directory for intermediate files`)
	sortCmd.Flags().IntP("max-open-files", "M", 400, `max number of open files`)
	sortCmd.Flags().BoolP("keep-tmp-dir", "k", false, `keep tmp dir`)
	sortCmd.Flags().BoolP("force", "", false, "overwrite tmp dir")
}
