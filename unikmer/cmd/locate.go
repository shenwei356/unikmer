// Copyright © 2018 Wei Shen <shenwei356@gmail.com>
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

	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/bio/seqio/fastx"
	"github.com/shenwei356/unikmer"
	"github.com/spf13/cobra"
)

// locateCmd represents
var locateCmd = &cobra.Command{
	Use:   "locate",
	Short: "locate Kmers in genome",
	Long: `locate Kmers in genome

Attention:
	1. output location is 1-based

`,
	Run: func(cmd *cobra.Command, args []string) {
		opt := getOptions(cmd)
		runtime.GOMAXPROCS(opt.NumCPUs)
		seq.ValidateSeq = false

		var err error

		var files []string
		infileList := getFlagString(cmd, "infile-list")
		if infileList != "" {
			files, err = getListFromFile(infileList)
			checkError(err)
		} else {
			files = getFileList(args)
		}

		checkFiles(files)
		if len(files) == 1 && isStdin(files[0]) {
			checkError(fmt.Errorf("stdin not supported, please give me .unik files"))
		}

		outFile := getFlagString(cmd, "out-prefix")
		circular := getFlagBool(cmd, "circular")

		genomeFile := getFlagString(cmd, "genome")

		// -----------------------------------------------------------------------

		var k int = -1
		var canonical bool

		var infh *bufio.Reader
		var r *os.File
		var reader *unikmer.Reader
		var nfiles = len(files)
		for i, file := range files {
			if isStdin(file) {
				log.Warningf("ignore stdin")
			}
			if opt.Verbose {
				log.Infof("pre-read file (%d/%d): %s", i+1, nfiles, file)
			}
			func() {
				infh, r, _, err = inStream(file)
				checkError(err)
				defer r.Close()

				reader, err = unikmer.NewReader(infh)
				checkError(err)

				if k == -1 {
					k = reader.K
					canonical = reader.Flag&unikmer.UNIK_CANONICAL > 0
					if opt.Verbose {
						if canonical {
							log.Infof("flag of canonical is on")
						} else {
							log.Infof("flag of canonical is off")
						}
					}
				} else if k != reader.K {
					checkError(fmt.Errorf("K (%d) of binary file '%s' not equal to previous K (%d)", reader.K, file, k))
				} else if (reader.Flag&unikmer.UNIK_CANONICAL > 0) != canonical {
					checkError(fmt.Errorf(`'canonical' flags not consistent, please check with "unikmer stats"`))
				}
			}()
		}

		// -----------------------------------------------------------------------

		m := make(map[uint64][]int, mapInitSize)

		var sequence, kmer, preKmer []byte
		var originalLen, l, end, e int
		var record *fastx.Record
		var fastxReader *fastx.Reader
		var kcode, preKcode unikmer.KmerCode
		var first bool
		var i, j, iters int
		var ok bool
		if opt.Verbose {
			log.Infof("read genome file: %s", genomeFile)
		}
		fastxReader, err = fastx.NewDefaultReader(genomeFile)
		checkError(err)
		for {
			record, err = fastxReader.Read()
			if err != nil {
				if err == io.EOF {
					break
				}
				checkError(err)
				break
			}

			if canonical {
				iters = 1
			} else {
				iters = 2
			}

			for j = 0; j < iters; j++ {
				if j == 0 { // sequence
					sequence = record.Seq.Seq

					if opt.Verbose {
						log.Infof("process sequence: %s", record.ID)
					}
				} else { // reverse complement sequence
					sequence = record.Seq.RevComInplace().Seq

					if opt.Verbose {
						log.Infof("process reverse complement sequence: %s", record.ID)
					}
				}

				originalLen = len(record.Seq.Seq)
				l = len(sequence)

				end = l - 1
				if end < 0 {
					end = 0
				}
				first = true
				for i = 0; i <= end; i++ {
					e = i + k
					if e > originalLen {
						if circular {
							e = e - originalLen
							kmer = sequence[i:]
							kmer = append(kmer, sequence[0:e]...)
						} else {
							break
						}
					} else {
						kmer = sequence[i : i+k]
					}

					if first {
						kcode, err = unikmer.NewKmerCode(kmer)
						first = false
					} else {
						kcode, err = unikmer.NewKmerCodeMustFromFormerOne(kmer, preKmer, preKcode)
					}
					if err != nil {
						checkError(fmt.Errorf("encoding '%s': %s", kmer, err))
					}
					preKmer, preKcode = kmer, kcode

					if canonical {
						kcode = kcode.Canonical()
					}

					if _, ok = m[kcode.Code]; !ok {
						m[kcode.Code] = make([]int, 0, 1)
					}
					if iters == 0 {
						m[kcode.Code] = append(m[kcode.Code], i)
					} else {
						m[kcode.Code] = append(m[kcode.Code], l-i-k)
					}
				}
			}
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

		var locs []int
		var loc int
		for i, file := range files {
			if isStdin(file) {
				log.Warningf("ignore stdin")
			}
			if opt.Verbose {
				log.Infof("process file (%d/%d): %s", i+1, nfiles, file)
			}
			func() {
				infh, r, _, err = inStream(file)
				checkError(err)
				defer r.Close()

				reader, err = unikmer.NewReader(infh)
				checkError(err)

				for {
					kcode, err = reader.Read()
					if err != nil {
						if err == io.EOF {
							break
						}
						checkError(err)
					}

					if locs, ok = m[kcode.Code]; ok {
						outfh.WriteString(kcode.String() + "\t")
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
	locateCmd.Flags().BoolP("circular", "", false, "circular genome")
	locateCmd.Flags().StringP("genome", "g", "", "genome in (gzipped) fasta file")
}
