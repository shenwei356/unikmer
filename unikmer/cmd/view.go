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
	"strings"

	"github.com/pkg/errors"
	"github.com/shenwei356/unikmer"
	"github.com/spf13/cobra"
)

// viewCmd represents
var viewCmd = &cobra.Command{
	Use:   "view",
	Short: "Read and output binary format to plain text",
	Long: `Read and output binary format to plain text

Attentions:
  1. The 'canonical/scaled/hashed' flags of all files should be consistent.
  2. Input files should ALL have or don't have taxid information.
  
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

		checkFileSuffix(opt, extDataFile, files...)

		outFile := getFlagString(cmd, "out-file")
		showCode := getFlagBool(cmd, "show-code")
		outFasta := getFlagBool(cmd, "fasta")
		outFastq := getFlagBool(cmd, "fastq")
		showCodeOnly := getFlagBool(cmd, "show-code-only")
		showTaxidOnly := getFlagBool(cmd, "show-taxid-only")
		genomes := getFlagStringSlice(cmd, "genome")
		providingGenomes := len(genomes) != 0

		showTaxid := getFlagBool(cmd, "show-taxid")
		if opt.IgnoreTaxid {
			showTaxid = false
		}

		outfh, gw, w, err := outStream(outFile, strings.HasSuffix(strings.ToLower(outFile), ".gz"), opt.CompressionLevel)
		checkError(err)
		defer func() {
			outfh.Flush()
			if gw != nil {
				gw.Close()
			}
			w.Close()
		}()

		var infh *bufio.Reader
		var r *os.File
		var reader0 *unikmer.Reader
		var canonical bool
		var hashed bool

		var k int = -1
		var hasTaxid bool

		var quality string

		var sequences [][]byte
		var hash2loc map[uint64][2]int // hash -> [seq idx, seq loc]
		var loc [2]int
		var ok bool

		var code uint64
		var kmer []byte
		var taxid uint32

		for _, file := range files {
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
					hasTaxid = !opt.IgnoreTaxid && reader.HasTaxidInfo()
					if showTaxid && !reader.HasTaxidInfo() {
						log.Warningf("flag -t/--show-taxid ignored when no taxids found in input")
					}
					if opt.IgnoreTaxid || !reader.HasTaxidInfo() {
						showTaxid = false
					}
					if providingGenomes {
						if !hashed {
							log.Warningf("-g/--genome ignored since k-mers not hashed")
						} else if !canonical {
							log.Warningf("-g/--genome ignored since 'canonical' flag is off")
						} else {
							if opt.Verbose {
								log.Infof("loading hash-k-mers mapping from genomes")
							}
							sequences, hash2loc, err = loadHash2Loc(genomes, k)
							checkError(errors.Wrap(err, "load hash-kmer mapping"))
							if opt.Verbose {
								log.Infof("%d hash-k-mers pairs from %d genomes loaded", len(hash2loc), len(sequences))
							}
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

				if outFastq {
					quality = strings.Repeat("g", reader.K)
				}

				for {
					code, taxid, err = reader.ReadCodeWithTaxid()
					if err != nil {
						if err == io.EOF {
							break
						}
						checkError(errors.Wrap(err, file))
					}

					if !hashed {
						kmer = unikmer.MustDecode(code, k)
					} else {
						if providingGenomes {
							if loc, ok = hash2loc[code]; ok {
								kmer = sequences[loc[0]][loc[1] : loc[1]+k]
							} else {
								kmer = []byte(fmt.Sprintf("%d", code))
								log.Warningf("fail to decode hash: %d, which is not found in given genomes", code)
							}
						} else {
							kmer = []byte(fmt.Sprintf("%d", code))
						}
					}

					if outFasta {
						if showTaxid {
							outfh.WriteString(fmt.Sprintf(">%d %d\n%s\n", code, taxid, kmer))
						} else {
							outfh.WriteString(fmt.Sprintf(">%d\n%s\n", code, kmer))
						}
					} else if outFastq {
						if showTaxid {
							outfh.WriteString(fmt.Sprintf("@%d %d\n%s\n+\n%s\n", code, taxid, kmer, quality))
						} else {
							outfh.WriteString(fmt.Sprintf("@%d\n%s\n+\n%s\n", code, kmer, quality))
						}
					} else if showTaxid {
						outfh.WriteString(fmt.Sprintf("%s\t%d\n", kmer, taxid))
					} else if showTaxidOnly {
						outfh.WriteString(fmt.Sprintf("%d\n", taxid))
					} else if showCodeOnly {
						outfh.WriteString(fmt.Sprintf("%d\n", code))
					} else if showCode {
						outfh.WriteString(fmt.Sprintf("%s\t%d\n", kmer, code))
					} else {
						outfh.WriteString(string(kmer) + "\n")
					}
				}

			}()
		}
	},
}

func init() {
	RootCmd.AddCommand(viewCmd)

	viewCmd.Flags().StringP("out-file", "o", "-", `out file ("-" for stdout, suffix .gz for gzipped out)`)
	viewCmd.Flags().BoolP("show-code", "n", false, `show encoded integer along with k-mer`)
	viewCmd.Flags().BoolP("show-code-only", "N", false, `only show encoded integers, faster than cutting from result of -n/--show-cde`)
	viewCmd.Flags().BoolP("fasta", "a", false, `output in FASTA format, with encoded integer as FASTA header`)
	viewCmd.Flags().BoolP("fastq", "q", false, `output in FASTQ format, with encoded integer as FASTQ header`)
	viewCmd.Flags().BoolP("show-taxid", "t", false, "show taxid")
	viewCmd.Flags().BoolP("show-taxid-only", "T", false, "show taxid only")
	viewCmd.Flags().StringSliceP("genome", "g", []string{}, "genomes in (gzipped) fasta file(s) for decoding hashed k-mers")
}
