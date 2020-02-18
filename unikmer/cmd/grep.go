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
	"strconv"
	"sync"

	"github.com/shenwei356/breader"
	"github.com/shenwei356/unikmer"
	"github.com/shenwei356/util/pathutil"
	"github.com/spf13/cobra"
)

// grepCmd represents
var grepCmd = &cobra.Command{
	Use:   "grep",
	Short: "Search k-mers from binary files",
	Long: `Search k-mers from binary files

Attentions:
  1. Canonical k-mers are used and outputed.
  2. Input files should ALL have or don't have taxid information.

Tips:
  1. Increase value of '-j' for better performance when dealing with
     lots of files, especially on SDD.

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
		queries := getFlagStringSlice(cmd, "query")
		queryFiles := getFlagStringSlice(cmd, "query-file")
		queryUnikFiles := getFlagStringSlice(cmd, "query-unik-file")
		queryWithTaxids := getFlagBool(cmd, "query-is-taxid")

		invertMatch := getFlagBool(cmd, "invert-match")
		degenerate := getFlagBool(cmd, "degenerate")

		mOutputs := getFlagBool(cmd, "multiple-outfiles")
		outdir := getFlagString(cmd, "out-dir")
		outSuffix := getFlagString(cmd, "out-suffix")
		force := getFlagBool(cmd, "force")

		sortKmers := getFlagBool(cmd, "sort")
		unique := getFlagBool(cmd, "unique")
		repeated := getFlagBool(cmd, "repeated")

		if (unique || repeated) && !sortKmers {
			log.Infof("flag -s/--sort is switched on when given -u/--unique or -d/--repeated")
			sortKmers = true
		}

		if len(queries) == 0 && len(queryFiles) == 0 && len(queryUnikFiles) == 0 {
			checkError(fmt.Errorf("one of flags -q/--query, -f/--query-file and -F/--query-unik-file needed"))
		}

		if mOutputs && !isStdin(outFile) {
			log.Warningf("flag -o/--out-prefix ignored when given -m/--multiple-outfiles")
		}

		var m map[uint64]struct{}
		var mt map[uint32]struct{}
		if queryWithTaxids {
			mt = make(map[uint32]struct{}, mapInitSize)
		} else {
			m = make(map[uint64]struct{}, mapInitSize)
		}

		k := -1

		// load k-mers from cli
		queryList := make([]string, 0, mapInitSize) // for plain k-mer text from -q and -f.
		for _, query := range queries {
			if query == "" {
				continue
			}
			if !queryWithTaxids {
				if k == -1 {
					k = len(query)
				} else if len(query) != k {
					checkError(fmt.Errorf("length of query sequence are inconsistent: (%d) != (%d): %s", len(query), k, query))
				}
			}
			queryList = append(queryList, query)
		}

		// load k-mers from file
		var nfiles int
		if len(queryFiles) != 0 {
			nfiles = len(queryFiles)
			var brdr *breader.BufferedReader
			var data interface{}
			var query string
			for i, queryFile := range queryFiles {
				if opt.Verbose {
					log.Infof("loading queries from k-mer file [%d/%d]: %s", i+1, nfiles, queryFile)
				}
				brdr, err = breader.NewDefaultBufferedReader(queryFile)
				checkError(err)
				for chunk := range brdr.Ch {
					checkError(chunk.Err)
					for _, data = range chunk.Data {
						query = data.(string)
						if !queryWithTaxids {
							if k == -1 {
								k = len(query)
							} else if len(query) != k {
								checkError(fmt.Errorf("length of query sequence are inconsistent: (%d) != (%d): %s", len(query), k, query))
							}
						}
						queryList = append(queryList, query)
					}
				}
			}
		}

		// encode k-mers or parse taxids
		var kcode unikmer.KmerCode
		var mer []byte
		var _queries [][]byte
		var q []byte
		var val uint64
		for _, query := range queryList {
			if queryWithTaxids {
				val, err = strconv.ParseUint(query, 10, 32)
				if err != nil {
					checkError(fmt.Errorf("query taxid should be positive integer in range of [1, %d]: %s", maxUint32, query))
				}

				mt[uint32(val)] = struct{}{}
				continue
			}
			if degenerate {
				_queries, err = extendDegenerateSeq([]byte(query))
				if err != nil {
					checkError(fmt.Errorf("fail to extend degenerate sequence '%s': %s", query, err))
				}
			} else {
				_queries = [][]byte{[]byte(query)}
			}

			for _, q = range _queries {
				kcode, err = unikmer.NewKmerCode(q)
				if err != nil {
					checkError(fmt.Errorf("fail to encode query '%s': %s", mer, err))
				}
				m[kcode.Canonical().Code] = struct{}{}
			}
		}

		var infh *bufio.Reader
		var r *os.File
		var reader *unikmer.Reader
		var flag int
		var canonical bool
		var taxid uint32

		// load k-mers/taxids from .unik files
		if len(queryUnikFiles) != 0 {
			nfiles = len(queryUnikFiles)
			for i, file := range queryUnikFiles {
				if opt.Verbose {
					log.Infof("loading queries from .unik file [%d/%d]: %s", i+1, nfiles, file)
				}

				flag = func() int {
					infh, r, _, err = inStream(file)
					checkError(err)
					defer r.Close()

					reader, err = unikmer.NewReader(infh)
					checkError(err)

					canonical = reader.IsCanonical()

					if queryWithTaxids && !reader.HasTaxidInfo() {
						checkError(fmt.Errorf("no taxids found in file: %s", file))
					}

					if k == -1 {
						k = reader.K
					} else if k != reader.K {
						checkError(fmt.Errorf("K (%d) of binary file '%s' not equal to previous K (%d)", reader.K, file, k))
					}

					for {
						kcode, taxid, err = reader.ReadWithTaxid()
						if err != nil {
							if err == io.EOF {
								break
							}
							checkError(err)
						}

						if queryWithTaxids {
							mt[taxid] = struct{}{}
							continue
						}
						if !canonical {
							kcode = kcode.Canonical()
						}
						m[kcode.Code] = struct{}{}
					}

					return flagContinue
				}()

				if flag == flagReturn {
					return
				} else if flag == flagBreak {
					break
				}
			}
		}

		if opt.Verbose {
			if queryWithTaxids {
				if len(mt) == 0 {
					log.Warningf("%d taxids loaded", len(mt))
					return
				}
				log.Infof("%d taxids loaded", len(mt))
			} else {
				if len(m) == 0 {
					log.Warningf("%d k-mers loaded", len(m))
					return
				}
				log.Infof("%d k-mers loaded", len(m))
			}
			log.Info()
		}

		// for faster query when only one taxid given
		var singleTaxidQuery, singleCodeQuery bool
		var theOneTaxid uint32
		var theOneCode uint64

		if queryWithTaxids {
			singleTaxidQuery = len(mt) == 1
			if singleTaxidQuery {
				for ot := range mt {
					theOneTaxid = ot
					break
				}
			}
		} else {
			singleCodeQuery = len(mt) == 1
			if singleCodeQuery {
				for oc := range m {
					theOneCode = oc
					break
				}
			}

		}

		////////////////////////////////////////////////////////////////////////////////

		var outfh *bufio.Writer
		var gw io.WriteCloser
		var w *os.File
		var writer *unikmer.Writer
		var hasTaxid bool

		if !mOutputs {
			// set global writer later
		} else {
			if outdir == "" {
				checkError(fmt.Errorf("out dir (flag -O/--out-dir) should not be empty"))
			}
			for _, file := range files {
				if isStdin(file) {
					checkError(fmt.Errorf("stdin detected, should not use -m/--mutliple-outfiles"))
				}
			}

			pwd, _ := os.Getwd()
			if outdir != "./" && outdir != "." && pwd != filepath.Clean(outdir) {
				existed, err := pathutil.DirExists(outdir)
				checkError(err)
				if existed {
					empty, err := pathutil.IsEmpty(outdir)
					checkError(err)
					if !empty {
						if force {
							checkError(os.RemoveAll(outdir))
						} else {
							log.Warningf("outdir not empty: %s, you can use --force to overwrite", outdir)
						}
					}
				} else {
					checkError(os.MkdirAll(outdir, 0777))
				}
			}
		}

		// -----------------------------------------------------------------------

		threads := opt.NumCPUs

		if threads > len(files) {
			threads = len(files)
		}
		if threads < 1 {
			threads = 1
		}

		if opt.Verbose {
			log.Infof("%d workers in position", threads)
		}

		var wg sync.WaitGroup
		tokens := make(chan int, threads)

		var codes []uint64
		var codesTaxids []unikmer.CodeTaxid
		if sortKmers {
			codes = make([]uint64, 0, mapInitSize)
			codesTaxids = make([]unikmer.CodeTaxid, 0, mapInitSize)
		}

		// read k-mers from goroutines
		var ns int
		var done chan int
		var chCodes chan uint64
		var chCodesTaxids chan unikmer.CodeTaxid
		var once sync.Once

		if !mOutputs {
			done = make(chan int)
			chCodes = make(chan uint64, threads)
			chCodesTaxids = make(chan unikmer.CodeTaxid, threads)
		}

		nfiles = len(files)
		for i, file := range files {
			tokens <- 1
			wg.Add(1)

			go func(i int, file string) {
				defer func() {
					<-tokens
					wg.Done()
				}()

				var infh *bufio.Reader
				var r *os.File
				var reader *unikmer.Reader

				var n int
				var _canonical bool
				var _hasGlobalTaxid bool
				var _isIncludeTaxid bool
				var _mustSort bool
				var _sorted bool
				var ok, hit bool

				if opt.Verbose {
					log.Infof("[file %d/%d] processing: %s", i+1, nfiles, file)
				}

				infh, r, _, err = inStream(file)
				checkError(err)
				defer r.Close()

				reader, err = unikmer.NewReader(infh)
				checkError(err)

				if !queryWithTaxids && k != reader.K {
					checkError(fmt.Errorf("K (%d) of binary file '%s' not equal to query K (%d)", reader.K, file, k))
				}

				_canonical = reader.IsCanonical()
				_hasGlobalTaxid = reader.HasGlobalTaxid()
				_isIncludeTaxid = reader.IsIncludeTaxid()
				_sorted = reader.IsSorted()

				// if the input files is already sorted, we don't have to sort again in mOutput mode.
				_mustSort = !reader.IsIncludeTaxid()

				if !mOutputs { // set global writer
					once.Do(func() {
						hasTaxid = !opt.IgnoreTaxid && reader.HasTaxidInfo()

						if !isStdout(outFile) {
							outFile += extDataFile
						}
						// the global writer
						outfh, gw, w, err = outStream(outFile, opt.Compress, opt.CompressionLevel)
						checkError(err)

						var mode uint32

						mode |= unikmer.UNIK_CANONICAL // forcing using canonical
						if sortKmers {
							mode |= unikmer.UNIK_SORTED
						} else if len(files) == 1 && reader.IsSorted() {
							// if the only input file is already sorted, we don't have to sort again.
							mode |= unikmer.UNIK_SORTED
						} else if opt.Compact {
							mode |= unikmer.UNIK_COMPACT
						}
						if hasTaxid {
							mode |= unikmer.UNIK_INCLUDETAXID
						}
						writer, err = unikmer.NewWriter(outfh, reader.K, mode)
						checkError(err)
						writer.SetMaxTaxid(maxUint32N(reader.GetTaxidBytesLength())) // follow reader

						go func() {
							if hasTaxid {
								for codeT := range chCodesTaxids {
									if sortKmers {
										codesTaxids = append(codesTaxids, codeT)
									} else {
										writer.WriteCodeWithTaxid(codeT.Code, codeT.Taxid)
										ns++
									}
								}
							} else {
								for code := range chCodes {
									if sortKmers {
										codes = append(codes, code)
									} else {
										writer.WriteCode(code)
										ns++
									}
								}
							}
							done <- 1
						}()
					})

					if !opt.IgnoreTaxid && reader.HasTaxidInfo() != hasTaxid {
						if reader.HasTaxidInfo() {
							checkError(fmt.Errorf(`taxid information not found in previous files, but found in this: %s`, file))
						} else {
							checkError(fmt.Errorf(`taxid information found in previous files, but missing in this: %s`, file))
						}
					}
				}

				var _writer *unikmer.Writer
				var _codes []uint64
				var _codesTaxids []unikmer.CodeTaxid
				var _outFile string

				if mOutputs {
					// write to it's own output file
					_outFile = filepath.Join(outdir, filepath.Base(file)+outSuffix+extDataFile)
					_outfh, _gw, _w, _err := outStream(_outFile, opt.Compress, opt.CompressionLevel)
					checkError(_err)
					defer func() {
						_outfh.Flush()
						if _gw != nil {
							_gw.Close()
						}
						_w.Close()
					}()

					var mode uint32
					mode |= unikmer.UNIK_CANONICAL
					if sortKmers {
						mode |= unikmer.UNIK_SORTED
					} else if reader.IsSorted() {
						mode |= unikmer.UNIK_SORTED
					} else if opt.Compact {
						mode |= unikmer.UNIK_COMPACT
					}
					if _isIncludeTaxid {
						mode |= unikmer.UNIK_INCLUDETAXID
					}
					_writer, err = unikmer.NewWriter(_outfh, reader.K, mode)
					checkError(err)
					_writer.SetMaxTaxid(maxUint32N(reader.GetTaxidBytesLength())) // follow reader
					if _hasGlobalTaxid {
						checkError(_writer.SetGlobalTaxid(reader.GetGlobalTaxid()))
					}

					if sortKmers && _mustSort {
						if _isIncludeTaxid {
							_codesTaxids = make([]unikmer.CodeTaxid, 0, mapInitSize)
						} else if _hasGlobalTaxid {
							_codes = make([]uint64, 0, mapInitSize)
						} else {
							_codes = make([]uint64, 0, mapInitSize)
						}
					}
					checkError(_writer.Flush())
				}

				var kcode unikmer.KmerCode
				var taxid uint32
				for {
					kcode, taxid, err = reader.ReadWithTaxid()
					if err != nil {
						if err == io.EOF {
							break
						}
						checkError(err)
					}

					if queryWithTaxids {
						if singleTaxidQuery {
							ok = taxid == theOneTaxid
						} else {
							_, ok = mt[taxid]
						}
					} else {
						if singleCodeQuery {
							ok = kcode.Code == theOneCode
							if _sorted && kcode.Code > theOneCode { // no need compare later codes
								break
							}
						} else {
							if !_canonical {
								kcode = kcode.Canonical()
							}
							_, ok = m[kcode.Code]
						}
					}

					if !invertMatch {
						hit = ok
					} else {
						hit = !ok
					}

					if !hit {
						continue
					}

					if mOutputs {
						if sortKmers && _mustSort {
							if _isIncludeTaxid {
								_codesTaxids = append(_codesTaxids, unikmer.CodeTaxid{Code: kcode.Code, Taxid: taxid})
							} else {
								_codes = append(_codes, kcode.Code)
							}
						} else {
							_writer.WriteCodeWithTaxid(kcode.Code, taxid)
							n++
						}
					} else {
						if hasTaxid {
							chCodesTaxids <- unikmer.CodeTaxid{Code: kcode.Code, Taxid: taxid}
						} else {
							chCodes <- kcode.Code
						}
					}
				}

				if !mOutputs {
					return
				}

				if sortKmers && _mustSort {
					if _isIncludeTaxid {
						if opt.Verbose {
							log.Infof("[file %d/%d] sorting %d k-mers", i+1, nfiles, len(_codesTaxids))
						}
						sort.Sort(unikmer.CodeTaxidSlice(_codesTaxids))
					} else {
						if opt.Verbose {
							log.Infof("[file %d/%d] sorting %d k-mers", i+1, nfiles, len(_codes))
						}
						sort.Sort(unikmer.CodeSlice(_codes))
					}

					if opt.Verbose {
						log.Infof("[file %d/%d] done sorting", i+1, nfiles)
					}

					if _isIncludeTaxid {
						if unique {
							var last uint64 = ^uint64(0)
							for _, codeT := range _codesTaxids {
								if codeT.Code == last {
									continue
								}
								last = codeT.Code
								n++
								_writer.WriteCodeWithTaxid(codeT.Code, codeT.Taxid)
							}
						} else if repeated {
							var last uint64 = ^uint64(0)
							var count int
							for _, codeT := range _codesTaxids {
								if codeT.Code == last {
									if count == 1 { // write once
										_writer.WriteCodeWithTaxid(codeT.Code, codeT.Taxid)
										n++
									}
									count++
								} else {
									last = codeT.Code
									count = 1
								}
							}
						} else {
							_writer.Number = int64(len(_codesTaxids))
							for _, codeT := range _codesTaxids {
								_writer.WriteCodeWithTaxid(codeT.Code, codeT.Taxid)
							}
							n = len(_codesTaxids)
						}
					} else {
						if unique {
							var last uint64 = ^uint64(0)
							for _, code := range _codes {
								if code == last {
									continue
								}
								last = code
								_writer.WriteCode(code)
								n++
							}
						} else if repeated {
							var last uint64 = ^uint64(0)
							var count int
							for _, code := range _codes {
								if code == last {
									if count == 1 { // write once
										_writer.WriteCode(code)
										n++
									}
									count++
								} else {
									last = code
									count = 1
								}
							}
						} else {
							_writer.Number = int64(len(_codes))
							for _, code := range _codes {
								_writer.WriteCode(code)
							}
							n = len(_codes)
						}
					}
				}

				checkError(_writer.Flush())
				if opt.Verbose {
					log.Infof("[file %d/%d] %d k-mers saved to %s", i+1, nfiles, n, _outFile)
				}

			}(i, file)
		}

		wg.Wait()

		if !mOutputs {
			close(chCodes)
			close(chCodesTaxids)
			<-done
		}

		if mOutputs {
			return
		}

		if sortKmers {
			if hasTaxid {
				if opt.Verbose {
					log.Infof("sorting %d k-mers", len(codesTaxids))
				}
				sort.Sort(unikmer.CodeTaxidSlice(codesTaxids))
			} else {
				if opt.Verbose {
					log.Infof("sorting %d k-mers", len(codes))
				}
				sort.Sort(unikmer.CodeSlice(codes))
			}
			if opt.Verbose {
				log.Infof("done sorting")
			}

			if hasTaxid {
				if unique {
					var last uint64 = ^uint64(0)
					for _, codeT := range codesTaxids {
						if codeT.Code == last {
							continue
						}
						last = codeT.Code
						writer.WriteCodeWithTaxid(codeT.Code, codeT.Taxid)
						ns++
					}
				} else if repeated {
					var last uint64 = ^uint64(0)
					var count int
					for _, codeT := range codesTaxids {
						if codeT.Code == last {
							if count == 1 { // write once
								writer.WriteCodeWithTaxid(codeT.Code, codeT.Taxid)
								ns++
							}
							count++
						} else {
							last = codeT.Code
							count = 1
						}
					}
				} else {
					writer.Number = int64(len(codesTaxids))
					for _, codeT := range codesTaxids {
						writer.WriteCodeWithTaxid(codeT.Code, codeT.Taxid)
					}
					ns = len(codesTaxids)
				}
			} else {
				if unique {
					var last uint64 = ^uint64(0)
					for _, code := range codes {
						if code == last {
							continue
						}
						last = code
						writer.WriteCode(code)
						ns++
					}
				} else if repeated {
					var last uint64 = ^uint64(0)
					var count int
					for _, code := range codes {
						if code == last {
							if count == 1 { // write once
								writer.WriteCode(code)
								ns++
							}
							count++
						} else {
							last = code
							count = 1
						}
					}
				} else {
					writer.Number = int64(len(codes))
					for _, code := range codes {
						writer.WriteCode(code)
					}
					ns = len(codes)
				}
			}
		}

		checkError(writer.Flush())

		checkError(outfh.Flush())
		if gw != nil {
			gw.Close()
		}
		w.Close()
		if opt.Verbose {
			log.Infof("%d k-mers saved to %s", ns, outFile)
		}
	},
}

func init() {
	RootCmd.AddCommand(grepCmd)

	grepCmd.Flags().StringP("out-prefix", "o", "-", `out file prefix ("-" for stdout)`)

	grepCmd.Flags().StringSliceP("query", "q", []string{""}, `query k-mers/taxids (multiple values delimted by comma supported)`)
	grepCmd.Flags().StringSliceP("query-file", "f", []string{""}, "query file (one k-mer/taxid per line)")
	grepCmd.Flags().StringSliceP("query-unik-file", "F", []string{""}, "query file in .unik format")
	grepCmd.Flags().BoolP("query-is-taxid", "t", false, "queries are taxids")

	grepCmd.Flags().BoolP("degenerate", "D", false, "query k-mers contains degenerate base")
	grepCmd.Flags().BoolP("invert-match", "v", false, "invert the sense of matching, to select non-matching records")

	grepCmd.Flags().BoolP("multiple-outfiles", "m", false, "write results into separated files for multiple input files")
	grepCmd.Flags().StringP("out-dir", "O", "unikmer-grep", "output directory")
	grepCmd.Flags().StringP("out-suffix", "S", grepDefaultOutSuffix, "output suffix")
	grepCmd.Flags().BoolP("force", "", false, "overwrite output directory")

	grepCmd.Flags().BoolP("sort", "s", false, helpSort)
	grepCmd.Flags().BoolP("unique", "u", false, `remove duplicated k-mers`)
	grepCmd.Flags().BoolP("repeated", "d", false, `only print duplicate k-mers`)

}

var grepDefaultOutSuffix = ".grep"
