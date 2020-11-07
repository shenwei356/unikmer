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
	"strings"
	"sync"

	"github.com/pkg/errors"
	"github.com/shenwei356/util/pathutil"

	"github.com/shenwei356/unikmer"
	"github.com/spf13/cobra"
)

var tsplitCmd = &cobra.Command{
	Use:   "tsplit",
	Short: "Split k-mers according to taxid",
	Long: `Split k-mers according to taxid

Attentions:
  1. The 'canonical' flags of all files should be consistent.
  2. Input files should ALL have taxid information.
  3. Input files should be sorted using 'unikmer sort'.
  4. All k-mers will loaded into RAM, for big input files,
     you can 'split' them first, 'tsplit' and then 'concat'
     for every taxid.
  
Tips:
  1. Increasing value of -j/--threads can accelerates splitting stage,
     in cost of more memory occupation.

`,
	Run: func(cmd *cobra.Command, args []string) {
		opt := getOptions(cmd)
		runtime.GOMAXPROCS(opt.NumCPUs)

		outdir := getFlagString(cmd, "out-dir")
		force := getFlagBool(cmd, "force")
		outPrefix := getFlagString(cmd, "out-prefix")

		if outPrefix == "" || strings.HasPrefix(outPrefix, ".") {
			checkError(fmt.Errorf(`-o/--out-prefix should not be empty or starting with "."`))
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

		var err error

		if outdir == "" {
			if isStdin(files[0]) {
				outdir = "stdin.tsplit"
			} else {
				outdir = files[0] + ".tsplit"
			}
		}
		pwd, _ := os.Getwd()
		if outdir != "./" && outdir != "." && pwd != filepath.Clean(outdir) {
			existed, err := pathutil.DirExists(outdir)
			checkError(errors.Wrap(err, outdir))
			if existed {
				empty, err := pathutil.IsEmpty(outdir)
				checkError(errors.Wrap(err, outdir))
				if !empty {
					if force {
						checkError(os.RemoveAll(outdir))
						checkError(os.MkdirAll(outdir, 0755))
					} else {
						log.Warningf("outdir not empty: %s, you can use --force to overwrite", outdir)
					}
				}
			} else {
				checkError(os.MkdirAll(outdir, 0755))
			}
		}

		m := make(map[uint32]*[]uint64, 1024) // taxid -> kmers

		var infh *bufio.Reader
		var r *os.File
		var reader0 *unikmer.Reader
		var codes *[]uint64
		var code uint64
		var taxid uint32
		var k int = -1
		var canonical bool
		var hashed bool
		var hasTaxid bool
		var mode uint32
		var flag int
		var nfiles = len(files)
		var ok bool
		var n int
		var maxTaxid uint32

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
					if !reader.IsSorted() {
						checkError(fmt.Errorf("input should be sorted: %s", file))
					}
					if canonical {
						mode |= unikmer.UnikCanonical
					}
					if hashed {
						mode |= unikmer.UnikHashed
					}
					mode |= unikmer.UnikSorted
					maxTaxid = maxUint32N(reader.GetTaxidBytesLength())
				} else {
					checkCompatibility(reader0, reader, file)
					if !opt.IgnoreTaxid && reader.HasTaxidInfo() != hasTaxid {
						if reader.HasTaxidInfo() {
							checkError(fmt.Errorf(`taxid information not found in previous files, but found in this: %s`, file))
						} else {
							checkError(fmt.Errorf(`taxid information found in previous files, but missing in this: %s`, file))
						}
					}
					if maxUint32N(reader.GetTaxidBytesLength()) > maxTaxid {
						maxTaxid = maxUint32N(reader.GetTaxidBytesLength())
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

					n++
					if codes, ok = m[taxid]; ok {
						*codes = append(*codes, code)
					} else {
						tmp := make([]uint64, 0, 1024)
						tmp = append(tmp, code)
						m[taxid] = &tmp
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

		if opt.Verbose {
			if n == 0 {
				log.Warningf("%d taxids loaded", n)
				return
			}
			log.Infof("%d taxids belonging to %d taxids loaded", n, len(m))
		}

		// ---------------------------------------
		threads := opt.NumCPUs
		if threads > len(m) {
			threads = len(m)
		}

		// just for counting total k-mers
		chN := make(chan int64, threads)
		var N int64 = 0
		done := make(chan int)
		go func() {
			tick := 0
			for n := range chN {
				N += n

				tick++
				if tick%100 == 0 {
					runtime.GC()
				}
			}
			done <- 1
		}()

		var wg sync.WaitGroup
		tokens := make(chan int, threads)

		var i int
		var ntaxids int = len(m)
		for _taxid, _codes := range m {
			wg.Add(1)
			tokens <- 1
			i++

			go func(taxid uint32, codes *[]uint64, i int) {
				defer func() {
					wg.Done()
					<-tokens
				}()

				_outFile := filepath.Join(outdir, fmt.Sprintf("%s.taxid-%d.k%d%s", outPrefix, taxid, k, extDataFile))
				_outfh, _gw, _w, _err := outStream(_outFile, opt.Compress, opt.CompressionLevel)
				checkError(_err)
				defer func() {
					_outfh.Flush()
					if _gw != nil {
						_gw.Close()
					}
					_w.Close()
				}()

				_writer, err := unikmer.NewWriter(_outfh, k, mode)
				checkError(errors.Wrap(err, _outFile))

				_writer.Number = int64(len(*codes))
				_writer.SetMaxTaxid(maxTaxid) // follow reader
				_writer.SetGlobalTaxid(taxid)

				for _, code := range *codes {
					_writer.WriteCode(code)
				}

				checkError(_writer.Flush())
				if opt.Verbose {
					log.Infof("[%d/%d] %d k-mers saved to %s", i, ntaxids, len(*codes), _outFile)
				}

				chN <- int64(len(*codes))

				*codes = make([]uint64, 1)
			}(_taxid, _codes, i)
		}

		// wait all k-mers being wrote to files
		wg.Wait()
		close(chN)
		<-done

		if opt.Verbose {
			log.Infof("%d taxids belonging to %d taxids saved to dir: %s", N, len(m), outdir)
		}
	},
}

func init() {
	RootCmd.AddCommand(tsplitCmd)

	tsplitCmd.Flags().StringP("out-prefix", "o", "tsplit", `out file prefix`)
	tsplitCmd.Flags().StringP("out-dir", "O", "", `output directory`)
	tsplitCmd.Flags().BoolP("force", "", false, `overwrite output directory`)
}
