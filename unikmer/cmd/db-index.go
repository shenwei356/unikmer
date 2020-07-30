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
	"math"
	"os"
	"path/filepath"
	"regexp"
	"runtime"
	"sort"
	"sync"

	"github.com/shenwei356/unikmer"
	"github.com/shenwei356/unikmer/index"
	"github.com/shenwei356/util/pathutil"
	"github.com/spf13/cobra"
)

// indexCmd represents
var indexCmd = &cobra.Command{
	Use:   "index",
	Short: "construct index from binary files",
	Long: `construct index from binary files

Attentions:
  0. All input files should be sorted, and output file is sorted.
  1. The 'canonical' flags of all files should be consistent.

`,
	Run: func(cmd *cobra.Command, args []string) {
		opt := getOptions(cmd)
		runtime.GOMAXPROCS(opt.NumCPUs)

		var err error

		// outFile := getFlagString(cmd, "out-prefix")

		nameRegexp := getFlagString(cmd, "name-regexp")
		extractName := nameRegexp != ""
		var reName *regexp.Regexp
		if extractName {
			if !regexp.MustCompile(`\(.+\)`).MatchString(nameRegexp) {
				checkError(fmt.Errorf(`value of -r (--name-regexp) must contains "(" and ")" to extract name`))
			}
			reName, err = regexp.Compile(nameRegexp)
			checkError(err)
		}

		sBlock := getFlagInt(cmd, "block-size")

		fpr := getFlagPositiveFloat64(cmd, "false-positive-rate")
		numHashes := getFlagPositiveInt(cmd, "num-hash")
		if numHashes > 255 {
			checkError(fmt.Errorf("value of -n/--num-hash too big: %d", numHashes))
		}

		outDir := getFlagString(cmd, "out-dir")
		force := getFlagBool(cmd, "force")

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

		existed, err := pathutil.DirExists(outDir)
		checkError(err)
		if existed {
			empty, err := pathutil.IsEmpty(outDir)
			checkError(err)
			if !empty {
				if force {
					checkError(os.RemoveAll(outDir))
				} else {
					checkError(fmt.Errorf("dir not empty: %s, choose another one or use --force to overwrite", outDir))
				}
			} else {
				checkError(os.RemoveAll(outDir))
			}
		}
		checkError(os.MkdirAll(outDir, 0777))

		// numSigs := CalcSignatureSize(61738843, nHash, fpr)
		// data := make([]byte, numSigs)

		// ------------------------------------------------------------------------------------
		// check unik files and read k-mers numbers
		fileInfos := make([]UnikFileInfo, 0, len(files))

		var infh *bufio.Reader
		var r *os.File
		var reader *unikmer.Reader
		// var code uint64
		var k int = -1
		var canonical bool
		var n int64
		// var nfiles = len(files)
		var name string
		var found [][]string
		names0 := make([]string, 0, len(files))
		for _, file := range files {
			// if opt.Verbose {
			// 	log.Infof("checking file (%d/%d): %s", i+1, nfiles, file)
			// }

			infh, r, _, err = inStream(file)
			checkError(err)

			reader, err = unikmer.NewReader(infh)
			checkError(err)

			if k == -1 {
				k = reader.K
				canonical = reader.IsCanonical()
			} else {
				if k != reader.K {
					checkError(fmt.Errorf("K (%d) of binary file '%s' not equal to previous K (%d)", reader.K, file, k))
				}
				if reader.IsCanonical() != canonical {
					checkError(fmt.Errorf(`'canonical' flags not consistent, please check with "unikmer stats"`))
				}
			}

			if reader.Number < 0 {
				checkError(fmt.Errorf("binary file not sorted or no k-mers number found: %s", file))
			}

			if extractName {
				found = reName.FindAllStringSubmatch(filepath.Base(file), -1)
				if len(found) > 0 {
					name = found[0][1]
				} else {
					name = filepath.Base(file)
				}
			} else {
				name = filepath.Base(file)
			}

			fileInfos = append(fileInfos, UnikFileInfo{Path: file, Name: name, Kmers: reader.Number})
			n += reader.Number
			r.Close()

			names0 = append(names0, name)
		}

		sort.Sort(UnikFileInfos(fileInfos))
		// fmt.Println(fileInfos)

		// ------------------------------------------------------------------------------------

		nFiles := len(fileInfos)
		if sBlock <= 0 {
			sBlock = int((math.Sqrt(float64(nFiles))+7)/8) * 8
		}
		if sBlock < 8 {
			sBlock = 8
		} else if sBlock > nFiles {
			sBlock = nFiles
		}

		if opt.Verbose {
			log.Infof("block size: %d", sBlock)
		}

		nIndexFiles := int((len(files) + sBlock - 1) / sBlock)
		indexFiles := make([]string, 0, nIndexFiles)
		ch := make(chan string, nIndexFiles)
		done := make(chan int)
		go func() {
			for f := range ch {
				indexFiles = append(indexFiles, f)
			}
			done <- 1
		}()

		var prefix string

		var b, j int
		var wg0 sync.WaitGroup
		tokens0 := make(chan int, opt.NumCPUs)
		tokensOpenFiles := make(chan int, 512)
		for i := 0; i < nFiles; i += sBlock {
			j = i + sBlock
			if j > nFiles {
				j = nFiles
			}

			b++

			prefix = fmt.Sprintf("[block #%03d]", b)
			if opt.Verbose {
				log.Infof("%s processing file #%d-#%d", prefix, i+1, j)
			}

			wg0.Add(1)
			tokens0 <- 1
			go func(files []UnikFileInfo, b int, prefix string) {
				var wg sync.WaitGroup
				tokens := make(chan int, opt.NumCPUs)

				var maxElements int64
				for _, info := range files {
					if maxElements < info.Kmers {
						maxElements = info.Kmers
					}
				}

				batchFiles := make([]string, 0, int((len(files)+7)/8))

				numSigs := CalcSignatureSize(uint64(maxElements), numHashes, fpr)
				if opt.Verbose {
					log.Infof("%s max number of k-mers: %d, number of signatures: %d", prefix, maxElements, numSigs)
				}

				var bb, jj int
				for ii := 0; ii < len(files); ii += 8 {
					jj = ii + 8
					if jj > len(files) {
						jj = len(files)
					}
					wg.Add(1)
					tokens <- 1
					bb++

					outFile := filepath.Join(outDir, fmt.Sprintf("block%03d_batch%03d%s", b, bb, extIndex))
					batchFiles = append(batchFiles, outFile)

					go func(_files []UnikFileInfo, bb int, maxElements int64, numSigs uint64, outFile string) {
						defer func() {
							wg.Done()
							<-tokens
						}()

						outfh, gw, w, err := outStream(outFile, true, opt.CompressionLevel)
						checkError(err)
						defer func() {
							outfh.Flush()
							if gw != nil {
								gw.Close()
							}
							w.Close()
						}()

						names := make([]string, 0, 8)
						for _, info := range _files {
							names = append(names, info.Name)
						}

						writer, err := index.NewWriter(outfh, k, canonical, uint8(numHashes), numSigs, names)
						checkError(err)
						defer func() {
							checkError(writer.Flush())
						}()

						sigs := make([]byte, numSigs)

						for _k, info := range _files {
							var infh *bufio.Reader
							var r *os.File
							var reader *unikmer.Reader
							var err error
							var code uint64

							tokensOpenFiles <- 1
							infh, r, _, err = inStream(info.Path)
							checkError(err)

							reader, err = unikmer.NewReader(infh)
							checkError(err)

							for {
								code, _, err = reader.ReadCodeWithTaxid()
								if err != nil {
									if err == io.EOF {
										break
									}
									checkError(err)
								}

								for _, loc := range hashLocations(code, numHashes, numSigs) {
									sigs[loc] |= 1 << (7 - _k)
								}

							}

							r.Close()

							<-tokensOpenFiles
						}

						checkError(writer.WriteBatch(sigs, len(sigs)))
						if opt.Verbose {
							log.Infof("%s batch #%03d: wrote %d signatures", prefix, bb, len(sigs))
						}

					}(files[ii:jj], bb, maxElements, numSigs, outFile)
				}

				wg.Wait()
				if opt.Verbose {
					log.Infof("%s merging %d index files", prefix, len(batchFiles))
				}

				blockFile := filepath.Join(outDir, fmt.Sprintf("block%03d%s", b, extIndex))
				checkError(MergeUnikIndex(opt, prefix, batchFiles, blockFile))

				if opt.Verbose {
					log.Infof("%s finished merging", prefix)
				}

				wg0.Done()
				ch <- filepath.Base(blockFile)
				<-tokens0
			}(fileInfos[i:j], b, prefix)
		}

		wg0.Wait()
		close(ch)
		<-done

		sort.Strings(indexFiles)
		dbInfo := NewUnikIndexDBInfo(int(index.Version), indexFiles)
		dbInfo.K = k
		dbInfo.Kmers = int(n)
		dbInfo.FPR = fpr
		dbInfo.Names = names0
		dbInfo.NumHashes = numHashes
		dbInfo.Canonical = canonical
		checkError(dbInfo.WriteTo(filepath.Join(outDir, dbInfoFile)))

		// ------------------------------------------------------------------------------------

		if opt.Verbose {
			log.Infof("index with %d k-mers saved to %s", n, outDir)
		}

	},
}

func init() {
	dbCmd.AddCommand(indexCmd)

	// indexCmd.Flags().StringP("out-prefix", "o", "-", `index file prefix ("-" for stdout)`)
	indexCmd.Flags().StringP("out-dir", "O", "", `output directory`)
	indexCmd.Flags().Float64P("false-positive-rate", "f", 0.3, `false positive rate of single bloom filter`)
	indexCmd.Flags().IntP("num-hash", "n", 1, `number of hashes`)
	indexCmd.Flags().IntP("block-size", "b", 0, `block size, default: sqrt(#.files)`)

	indexCmd.Flags().BoolP("force", "", false, "overwrite tmp dir")
	indexCmd.Flags().StringP("name-regexp", "r", "", "regular expression for extract name from file name")

}
