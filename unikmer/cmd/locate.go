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
	"strings"

	"github.com/pkg/errors"
	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/bio/seqio/fastx"
	"github.com/shenwei356/unikmer"
	"github.com/spf13/cobra"
)

var locateCmd = &cobra.Command{
	Use:   "locate",
	Short: "Locate k-mers in genome",
	Long: `Locate k-mers in genome

Attention:
  1. All files should have the 'canonical' flag.
  2. Output location is 1-based.

`,
	Run: func(cmd *cobra.Command, args []string) {
		opt := getOptions(cmd)
		runtime.GOMAXPROCS(opt.NumCPUs)
		seq.ValidateSeq = false

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

		if len(files) == 1 && isStdin(files[0]) {
			checkError(fmt.Errorf("stdin not supported, please give me .unik files"))
		}

		outFile := getFlagString(cmd, "out-prefix")

		genomeFile := getFlagNonEmptyString(cmd, "genome")

		showHash := getFlagBool(cmd, "show-hash")

		// -----------------------------------------------------------------------

		var k int = -1
		var canonical bool
		var hashed bool

		var infh *bufio.Reader
		var r *os.File
		var reader0 *unikmer.Reader
		var nfiles = len(files)
		for i, file := range files {
			if isStdin(file) {
				log.Warningf("ignoring stdin")
			}
			if opt.Verbose {
				log.Infof("pre-reading file (%d/%d): %s", i+1, nfiles, file)
			}
			func() {
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
					if !canonical {
						checkError(fmt.Errorf("%s: 'canonical' flag is needed", file))
					}
				} else {
					checkCompatibility(reader0, reader, file)
				}

			}()
		}

		// -----------------------------------------------------------------------

		m := make(map[uint64][]int, mapInitSize)
		var hash2loc map[uint64][2]int // hash -> [seq idx, seq loc]
		var sequences [][]byte
		if hashed && !showHash {
			sequences = make([][]byte, 0, 8)
			hash2loc = make(map[uint64][2]int, mapInitSize)
		}

		var fastxReader *fastx.Reader
		var record *fastx.Record
		var iter *unikmer.Iterator
		var code uint64
		var ok bool
		if opt.Verbose {
			log.Infof("reading genome file: %s", genomeFile)
		}
		fastxReader, err = fastx.NewDefaultReader(genomeFile)
		checkError(errors.Wrap(err, genomeFile))
		var seqIdx, kmerIdx int
		for {
			record, err = fastxReader.Read()
			if err != nil {
				if err == io.EOF {
					break
				}
				checkError(errors.Wrap(err, genomeFile))
				break
			}
			// using ntHash
			if hashed {
				iter, err = unikmer.NewHashIterator(record.Seq, k, true)
			} else {
				iter, err = unikmer.NewKmerIterator(record.Seq, k, true)
			}
			checkError(errors.Wrapf(err, "seq: %s", record.Name))

			kmerIdx = 0
			if hashed {
				if !showHash {
					sequences = append(sequences, record.Seq.Clone().Seq)
				}
				for {
					code, ok = iter.NextHash()
					if !ok {
						break
					}

					if _, ok = m[code]; !ok {
						m[code] = make([]int, 0, 1)

						if hashed && !showHash {
							hash2loc[code] = [2]int{seqIdx, kmerIdx}
						}
					}
					m[code] = append(m[code], kmerIdx)

					kmerIdx++
				}
			} else {
				for {
					code, ok, err = iter.NextKmer()
					checkError(errors.Wrapf(err, "seq: %s", record.Name))
					if !ok {
						break
					}

					if _, ok = m[code]; !ok {
						m[code] = make([]int, 0, 1)
					}
					m[code] = append(m[code], kmerIdx)

					kmerIdx++
				}
			}

			seqIdx++
		}
		if opt.Verbose {
			log.Infof("finished reading genome file: %s", genomeFile)
		}

		// -----------------------------------------------------------------------

		outfh, gw, w, err := outStream(outFile, strings.HasSuffix(strings.ToLower(outFile), ".gz"), opt.CompressionLevel)
		checkError(err)
		defer func() {
			outfh.Flush()
			if gw != nil {
				gw.Close()
			}
			w.Close()
		}()

		var reader *unikmer.Reader
		var locs []int
		var loc int
		var seqKmerIdx [2]int
		for i, file := range files {
			if isStdin(file) {
				log.Warningf("ignoring stdin")
			}
			if opt.Verbose {
				log.Infof("processing file (%d/%d): %s", i+1, nfiles, file)
			}
			func() {
				infh, r, _, err = inStream(file)
				checkError(err)
				defer r.Close()

				reader, err = unikmer.NewReader(infh)
				checkError(errors.Wrap(err, file))

				for {
					code, _, err = reader.ReadCodeWithTaxid()
					if err != nil {
						if err == io.EOF {
							break
						}
						checkError(errors.Wrap(err, file))
					}

					if locs, ok = m[code]; ok {
						// locs = uniqInts(locs)
						sort.Ints(locs)
						if hashed {
							if showHash {
								outfh.WriteString(fmt.Sprintf("%d\t", code))
							} else {
								seqKmerIdx = hash2loc[code]
								outfh.WriteString(string(sequences[seqKmerIdx[0]][seqKmerIdx[1]:seqKmerIdx[1]+k]) + "\t")
							}
						} else {
							outfh.WriteString(unikmer.KmerCode{Code: code, K: k}.String() + "\t")
						}
						outfh.WriteString(fmt.Sprintf("%d", locs[0]+1))
						for _, loc = range locs[1:] {
							outfh.WriteString(fmt.Sprintf(",%d", loc+1))
						}
						outfh.WriteString("\n")
					}
				}
			}()
		}
	},
}

func init() {
	RootCmd.AddCommand(locateCmd)

	locateCmd.Flags().StringP("out-prefix", "o", "-", `out file prefix ("-" for stdout)`)
	locateCmd.Flags().StringP("genome", "g", "", "genome in (gzipped) fasta file")
	locateCmd.Flags().BoolP("show-hash", "H", false, "show hash instead of kmer for hashed kmer")
}
