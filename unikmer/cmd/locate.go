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
  0. All files should have the 'canonical' flag.
  1. The 'canonical/scaled/hashed' flags of all files should be consistent.
  2. Output is BED6 format.
  3. When using experimental flag --circular, leading subsequence of k-1 bp
     is appending to end of sequence. End position of k-mers that crossing
     sequence end would be greater than sequence length.

`,
	Run: func(cmd *cobra.Command, args []string) {
		opt := getOptions(cmd)
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

		genomes := getFlagStringSlice(cmd, "genome")
		if len(genomes) == 0 {
			checkError(fmt.Errorf("flag -g/--genome needed"))
		}

		circular := getFlagBool(cmd, "circular")

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

		m := make(map[uint64][][2]int, mapInitSize) // code -> locs
		var sequences [][]byte
		var ids [][]byte

		sequences = make([][]byte, 0, 8)
		ids = make([][]byte, 0, 8)

		var fastxReader *fastx.Reader
		var record *fastx.Record
		var iter *unikmer.Iterator
		var code uint64
		var ok bool
		var seqIdx int

		for _, file := range genomes {
			if opt.Verbose {
				log.Infof("reading genome file: %s", file)
			}
			fastxReader, err = fastx.NewDefaultReader(file)
			checkError(errors.Wrap(err, file))

			for {
				record, err = fastxReader.Read()
				if err != nil {
					if err == io.EOF {
						break
					}
					checkError(errors.Wrap(err, file))
					break
				}
				// using ntHash
				if hashed {
					iter, err = unikmer.NewHashIterator(record.Seq, k, true, circular)
				} else {
					iter, err = unikmer.NewKmerIterator(record.Seq, k, true, circular)
				}
				if err != nil {
					if err == unikmer.ErrShortSeq {
						if opt.Verbose {
							log.Infof("ignore short seq: %s", record.Name)
						}
						continue
					} else {
						checkError(errors.Wrapf(err, "seq: %s", record.Name))
					}
				}

				seqClone := record.Seq.Clone().Seq
				if circular {
					seqClone = append(seqClone, seqClone[0:k-1]...)
				}
				sequences = append(sequences, seqClone)
				ids = append(ids, record.ID)

				for {
					code, ok, err = iter.Next()
					if !hashed && err != nil {
						checkError(errors.Wrapf(err, "%s: %s: %s", file, record.Name, sequences[iter.Index():iter.Index()+k]))
					}
					if !ok {
						break
					}

					if _, ok = m[code]; !ok {
						m[code] = make([][2]int, 0, 1)
					}
					m[code] = append(m[code], [2]int{seqIdx, iter.Index()})

				}

				seqIdx++
			}
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
		var locs [][2]int
		var loc [2]int
		var j int
		var kmer []byte
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
						for _, loc = range locs {
							i, j = loc[0], loc[1]

							kmer = sequences[i][j : j+k]

							outfh.WriteString(fmt.Sprintf("%s\t%d\t%d\t%s\t0\t.\n",
								ids[i], j, j+k, kmer))
						}
					}
				}
			}()
		}
	},
}

func init() {
	RootCmd.AddCommand(locateCmd)

	locateCmd.Flags().StringP("out-prefix", "o", "-", `out file prefix ("-" for stdout)`)
	locateCmd.Flags().StringSliceP("genome", "g", []string{}, "genomes in (gzipped) fasta file(s)")
	locateCmd.Flags().BoolP("circular", "", false, `circular genome. type "unikmer locate -h" for details`)
}
