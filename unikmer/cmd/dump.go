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
	"fmt"
	"strconv"
	"strings"

	"github.com/pkg/errors"
	"github.com/shenwei356/breader"
	"github.com/shenwei356/kmers"
	"github.com/shenwei356/unik/v5"

	"github.com/spf13/cobra"
	"github.com/will-rowe/nthash"
)

var dumpCmd = &cobra.Command{
	Use:   "dump",
	Short: "Convert plain k-mer text to binary format",
	Long: `Convert plain k-mer text to binary format


Attentions:
  1. Input should be one k-mer per line, or tab-delimited two columns
     with a k-mer and it's taxid.
  2. You can also assign a global taxid with flag -t/--taxid.

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

		outFile := getFlagString(cmd, "out-prefix")
		unique := getFlagBool(cmd, "unique")
		canonical := getFlagBool(cmd, "canonical")
		canonicalOnly := getFlagBool(cmd, "canonical-only")
		sortedKmers := getFlagBool(cmd, "sorted")
		taxid := getFlagUint32(cmd, "taxid")
		hashed := getFlagBool(cmd, "hash")          // to compue the hash values of k-mers
		hashedAlready := getFlagBool(cmd, "hashed") // what given are hash values

		if hashed && canonicalOnly {
			checkError(fmt.Errorf("flag -H/--hash and -k/--canonical-only are not compatible"))
		}

		var k int = -1
		if hashedAlready {
			canonical = true

			k = getFlagNonNegativeInt(cmd, "kmer-len")
			if k == 0 {
				checkError(fmt.Errorf("flag -k/--kmer-len should be given when --hashed given"))
			}
		}

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

		var m map[uint64]struct{}
		if unique {
			m = make(map[uint64]struct{}, mapInitSize)
		}

		var l int
		var reader *breader.BufferedReader
		var chunk breader.Chunk
		var data interface{}
		var line string
		var linebytes []byte
		var kcode, kcodeC kmers.KmerCode
		var ok bool
		var n int64

		var includeTaxid bool
		var items []string
		var tmp uint64
		var _taxid uint32
		var once bool = true
		hasGlobalTaxid := taxid > 0
		var hasher *nthash.NTHi
		var hash uint64

		for _, file := range files {
			reader, err = breader.NewDefaultBufferedReader(file)
			checkError(errors.Wrap(err, file))

			for chunk = range reader.Ch {
				checkError(chunk.Err)
				for _, data = range chunk.Data {
					line = data.(string)
					l = len(line)

					if l == 0 {
						continue
					}

					if !hashedAlready {
						if k == -1 {
							if strings.Index(line, "\t") > 0 {
								includeTaxid = true
								if hasGlobalTaxid {
									log.Warningf("since input has more than one column, global taxid is ignored")
								}
							} else if hasGlobalTaxid {
								_taxid = taxid
								k = l
							} else {
								k = l
							}
						} else if !includeTaxid && l != k {
							checkError(fmt.Errorf("K-mer length mismatch, previous: %d, current: %d. %s", k, l, line))
						}
					} else {
						if once {
							if strings.Index(line, "\t") > 0 {
								includeTaxid = true
								if hasGlobalTaxid {
									log.Warningf("since input has more than one column, global taxid is ignored")
								}
							} else if hasGlobalTaxid {
								_taxid = taxid
							}
							once = false
						}
					}

					if includeTaxid {
						items = strings.Split(line, "\t")
						if len(items) < 2 {
							checkError(fmt.Errorf("inconsistent two column tabular format"))
						}
						line = items[0]
						l = len(line)

						if !hashedAlready {
							if once {
								k = len(line)
								once = false
							} else {
								if l != k {
									checkError(fmt.Errorf("K-mer length mismatch, previous: %d, current: %d. %s", k, l, line))
								}
							}
						}

						tmp, err = strconv.ParseUint(items[1], 10, 32)
						if err != nil {
							checkError(fmt.Errorf("query taxid (2nd column) should be positive integer in range of [1, %d]: %s", maxUint32, items[1]))
						}

						_taxid = uint32(tmp)
					}

					if writer == nil {
						var mode uint32
						if sortedKmers {
							mode |= unik.UnikSorted
						} else if opt.Compact && !hashed {
							mode |= unik.UnikCompact
						}
						if canonical || canonicalOnly {
							mode |= unik.UnikCanonical
						}
						if includeTaxid {
							mode |= unik.UnikIncludeTaxID
						}
						if hashed || hashedAlready {
							mode |= unik.UnikHashed
						}
						writer, err = unik.NewWriter(outfh, k, mode)
						if err != nil {
							checkError(errors.Wrap(err, outFile))
						}
						writer.SetMaxTaxid(opt.MaxTaxid)
						if !includeTaxid && hasGlobalTaxid {
							checkError(writer.SetGlobalTaxid(taxid))
						}
					}

					if hashedAlready {
						hash, err = strconv.ParseUint(line, 10, 64)
						if err != nil {
							checkError(err)
						}

						if unique {
							if _, ok = m[hash]; !ok {
								m[hash] = struct{}{}
								checkError(writer.WriteCode(hash))
								if includeTaxid {
									checkError(writer.WriteTaxid(_taxid))
								}
								n++
							}
						} else {
							checkError(writer.WriteCode(hash))
							if includeTaxid {
								checkError(writer.WriteTaxid(_taxid))
							}
							n++
						}

						continue
					}

					if hashed {
						linebytes = []byte(line)
						hasher, err = nthash.NewHasher(&linebytes, uint(k))
						if err != nil {
							checkError(errors.Wrap(err, line))
						}
						// for hash = range hasher.Hash(canonical) {
						// 	break
						// }
						hash, _ = hasher.Next(canonical)

						if unique {
							if _, ok = m[hash]; !ok {
								m[hash] = struct{}{}
								checkError(writer.WriteCode(hash))
								if includeTaxid {
									checkError(writer.WriteTaxid(_taxid))
								}
								n++
							}
						} else {
							checkError(writer.WriteCode(hash))
							if includeTaxid {
								checkError(writer.WriteTaxid(_taxid))
							}
							n++
						}

						continue
					}

					kcode, err = kmers.NewKmerCode([]byte(line))
					if err != nil {
						checkError(fmt.Errorf("fail to encode '%s': %s", line, err))
					}

					if canonicalOnly {
						kcodeC = kcode.Canonical()
						if kcode.Code != kcodeC.Code {
							continue
						}
						kcode = kcodeC
					} else if canonical {
						kcode = kcode.Canonical()
					}

					if unique {
						if _, ok = m[kcode.Code]; !ok {
							m[kcode.Code] = struct{}{}
							checkError(writer.WriteCode(kcode.Code))
							if includeTaxid {
								checkError(writer.WriteTaxid(_taxid))
							}
							n++
						}
					} else {
						checkError(writer.WriteCode(kcode.Code))
						if includeTaxid {
							checkError(writer.WriteTaxid(_taxid))
						}
						n++
					}
				}
			}
		}

		checkError(writer.Flush())
		if opt.Verbose {
			log.Infof("%d unique k-mers saved to %s", n, outFile)
		}
	},
}

func init() {
	RootCmd.AddCommand(dumpCmd)

	dumpCmd.Flags().StringP("out-prefix", "o", "-", `out file prefix ("-" for stdout)`)
	dumpCmd.Flags().BoolP("unique", "u", false, `remove duplicate k-mers`)
	dumpCmd.Flags().BoolP("canonical", "K", false, "save the canonical k-mers")
	dumpCmd.Flags().BoolP("canonical-only", "O", false, "only save the canonical k-mers. This flag overides -K/--canonical")
	dumpCmd.Flags().BoolP("sorted", "s", false, "input k-mers are sorted")
	dumpCmd.Flags().Uint32P("taxid", "t", 0, "global taxid")
	dumpCmd.Flags().BoolP("hash", "H", false, `save hash of k-mer, automatically on for k>32. This flag overides global flag -c/--compact`)

	dumpCmd.Flags().BoolP("hashed", "", false, `giving hash values of k-mers, This flag overides global flag -c/--compact`)
	dumpCmd.Flags().IntP("kmer-len", "k", 0, "k-mer length")
}
