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

	"github.com/pkg/errors"
	"github.com/shenwei356/unikmer"
	"github.com/shenwei356/unikmer/index"
	"github.com/shenwei356/util/bytesize"
	"github.com/shenwei356/util/pathutil"
	"github.com/spf13/cobra"
)

var indexCmd = &cobra.Command{
	Use:   "index",
	Short: "Construct index from binary files",
	Long: `Construct index from binary files

Attentions:
  0. All input files should be sorted.
  1. The 'canonical' flags of all files should be consistent.
  2. Increase value of -j/--threads for acceleratation in cost of more
     memory occupation and I/O pressure, sqrt(#cpus) is recommended.
     #threads * #threads files are simultaneously opened, and max number
     of opened files is limited by flag -F/--max-open-files.
  3. Value of block size -b/--block-size better be multiple of 64.
  4. Use --dry-run to adjust parameters and see final #index files 
     and total file size.

`,
	Run: func(cmd *cobra.Command, args []string) {
		opt := getOptions(cmd)
		runtime.GOMAXPROCS(opt.NumCPUs)

		dryRun := getFlagBool(cmd, "dry-run")
		if dryRun {
			opt.Verbose = true
		}

		// block-max-kmers-t1
		kmerThreshold8Str := getFlagString(cmd, "block-max-kmers-t1")
		kmerThreshold8Float, err := bytesize.Parse([]byte(kmerThreshold8Str))
		if err != nil {
			checkError(fmt.Errorf("invalid size: %s", kmerThreshold8Float))
		}
		kmerThreshold8 := int64(kmerThreshold8Float)
		if kmerThreshold8 <= 0 {
			checkError(fmt.Errorf("value of flag -m/--block-max-kmers-t1 should be positive: %d", kmerThreshold8))
		}

		// block-max-kmers-t2
		kmerThresholdSStr := getFlagString(cmd, "block-max-kmers-t2")
		kmerThresholdSFloat, err := bytesize.Parse([]byte(kmerThresholdSStr))
		if err != nil {
			checkError(fmt.Errorf("invalid size: %s", kmerThresholdSFloat))
		}
		kmerThresholdS := int64(kmerThresholdSFloat)
		if kmerThresholdS <= 0 {
			checkError(fmt.Errorf("value of flag -M/--block-max-kmers-t2 should be positive: %d", kmerThresholdS))
		}

		if kmerThreshold8 >= kmerThresholdS {
			checkError(fmt.Errorf("value of flag -m/--block-max-kmers-t1 (%d) should be small than -M/--block-max-kmers-t2 (%d)", kmerThreshold8, kmerThresholdS))
		}

		// -------------------------------------------------------

		nameRegexp := getFlagString(cmd, "name-regexp")
		extractName := nameRegexp != ""
		var reName *regexp.Regexp
		if extractName {
			if !regexp.MustCompile(`\(.+\)`).MatchString(nameRegexp) {
				checkError(fmt.Errorf(`value of -r (--name-regexp) must contains "(" and ")" to extract name`))
			}
			reName, err = regexp.Compile(nameRegexp)
			checkError(errors.Wrapf(err, "parsing regular expression: %s", nameRegexp))
		}

		sBlock := getFlagInt(cmd, "block-size")

		fpr := getFlagPositiveFloat64(cmd, "false-positive-rate")
		numHashes := getFlagPositiveInt(cmd, "num-hash")
		if numHashes > 255 {
			checkError(fmt.Errorf("value of -n/--num-hash too big: %d", numHashes))
		}
		maxOpenFiles := getFlagPositiveInt(cmd, "max-open-files")

		outDir := getFlagString(cmd, "out-dir")
		if outDir == "" {
			checkError(fmt.Errorf("value of -O/--out-dir can not be empty"))
		}

		force := getFlagBool(cmd, "force")

		if opt.Verbose {
			log.Infof("number of CPUs to use: %d", opt.NumCPUs)
			log.Infof("number of hashes: %d", numHashes)
			log.Infof("false positive rate: %f", fpr)
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
		if len(files) < 2 {
			checkError(fmt.Errorf("at least two .unik files needed"))
		}

		existed, err := pathutil.DirExists(outDir)
		checkError(errors.Wrapf(err, "check output dir: %s", outDir))
		if existed {
			empty, err := pathutil.IsEmpty(outDir)
			checkError(errors.Wrapf(err, "check output dir: %s", outDir))
			if !empty {
				if force {
					if !dryRun {
						checkError(os.RemoveAll(outDir))
					}
				} else {
					if !dryRun {
						checkError(fmt.Errorf("dir not empty: %s, choose another one or use --force to overwrite", outDir))
					}
				}
			} else {
				if !dryRun {
					checkError(os.RemoveAll(outDir))
				}
			}
		}
		if !dryRun {
			checkError(os.MkdirAll(outDir, 0777))
		}

		// ------------------------------------------------------------------------------------
		// check unik files and read k-mers numbers

		if opt.Verbose {
			log.Info("checking .unik files ...")
		}

		fileInfos := make([]UnikFileInfo, 0, len(files))

		var k int = -1
		var canonical bool
		var n int64
		var nfiles = len(files)

		getInfo := func(file string, first bool) UnikFileInfo {
			infh, r, _, err := inStream(file)
			checkError(err)

			reader, err := unikmer.NewReader(infh)
			checkError(errors.Wrap(err, file))

			if first {
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
			var name string
			if extractName {
				found := reName.FindAllStringSubmatch(filepath.Base(file), -1)
				if len(found) > 0 {
					name = found[0][1]
				} else {
					name = filepath.Base(file)
				}
			} else {
				name = filepath.Base(file)
			}

			checkError(r.Close())
			return UnikFileInfo{Path: file, Name: name, Kmers: reader.Number}
		}

		// fisrt file
		if opt.Verbose {
			log.Infof("checking file: %d/%d", 1, nfiles)
		}
		file := files[0]
		info := getInfo(file, true)
		fileInfos = append(fileInfos, info)
		namesMap := make(map[string]interface{}, nfiles)
		namesMap[info.Name] = struct{}{}

		// left files
		var wgGetInfo sync.WaitGroup
		chInfos := make(chan UnikFileInfo, opt.NumCPUs)
		tokensGetInfo := make(chan int, opt.NumCPUs)
		doneGetInfo := make(chan int)
		go func() {
			var ok bool
			for info := range chInfos {
				fileInfos = append(fileInfos, info)
				n += info.Kmers

				if _, ok = namesMap[info.Name]; ok {
					log.Warningf("duplicated name: %s", info.Name)
				} else {
					namesMap[info.Name] = struct{}{}
				}
			}
			doneGetInfo <- 1
		}()
		for i, file := range files[1:] {
			if opt.Verbose {
				if i < 98 || (i+2)%100 == 0 {
					log.Infof("checking file: %d/%d", i+2, nfiles)
				}
			}
			wgGetInfo.Add(1)
			tokensGetInfo <- 1
			go func(file string) {
				defer func() {
					wgGetInfo.Done()
					<-tokensGetInfo
				}()
				chInfos <- getInfo(file, false)
			}(file)
		}

		wgGetInfo.Wait()
		close(chInfos)
		<-doneGetInfo

		if opt.Verbose {
			log.Infof("finished checking %d .unik files", nfiles)
		}

		sort.Sort(UnikFileInfos(fileInfos))

		names0 := make([]string, 0, len(files))
		sizes0 := make([]uint64, 0, len(files))
		for _, info := range fileInfos {
			n += info.Kmers
			names0 = append(names0, info.Name)
			sizes0 = append(sizes0, uint64(info.Kmers))
		}

		if dryRun {
			log.Infof("names:")
			for _, info := range fileInfos {
				log.Infof("name: %s, #k-mers: %d, file: %s, ", info.Name, info.Kmers, info.Path)
			}
		}

		// ------------------------------------------------------------------------------------
		// begin creating index

		if opt.Verbose {
			log.Infof("starting indexing ...")
		}
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
			log.Infof("block-max-kmers-threshold 1: %s", kmerThreshold8Str)
			log.Infof("block-max-kmers-threshold 2: %s", kmerThresholdSStr)
		}

		nIndexFiles := int((len(files) + sBlock - 1) / sBlock) // may be more if using -m and -M
		indexFiles := make([]string, 0, nIndexFiles)

		ch := make(chan string, nIndexFiles)
		done := make(chan int)
		go func() {
			for f := range ch {
				indexFiles = append(indexFiles, f)
			}
			done <- 1
		}()

		var fileSize float64
		chFileSize := make(chan float64, nIndexFiles)
		doneFileSize := make(chan int)
		go func() {
			for f := range chFileSize {
				fileSize += f
			}
			doneFileSize <- 1
		}()

		var prefix string

		var b int
		var wg0 sync.WaitGroup
		maxConc := opt.NumCPUs
		if dryRun {
			maxConc = 1 // just for loging in order
		}
		tokens0 := make(chan int, maxConc)
		tokensOpenFiles := make(chan int, maxOpenFiles)

		sBlock0 := sBlock

		batch := make([]*UnikFileInfo, 0, sBlock)
		var flag8, flag bool
		var lastInfo *UnikFileInfo
		for i := 0; i <= nFiles; i++ {
			if i == nFiles { // process lastInfo
				// process files in batch or the last file
				if flag {
					batch = append(batch, lastInfo)
				}
			} else {
				info := fileInfos[i]

				if flag || flag8 {
					if flag {
						batch = append(batch, lastInfo)
						lastInfo = &info
					} else if info.Kmers > kmerThresholdS {
						// meet a very big file the first time
						flag = true      // mark
						lastInfo = &info // leave this file process in the next round
						// and we have to process files aleady in batch
					} else { // flag8 && !flag
						if lastInfo != nil { // right after found the first file > kmerThreshold8
							batch = append(batch, lastInfo)
							lastInfo = nil
						}

						batch = append(batch, &info)
						if len(batch) < sBlock { // not filled
							continue
						}
					}
				} else if info.Kmers > kmerThreshold8 {
					if info.Kmers > kmerThresholdS {
						// meet a very big file the first time
						flag = true      // mark
						lastInfo = &info // leave this file process in the next round
						// and we have to process files aleady in batch
					} else {
						// meet a big file > kmerThreshold8
						sBlock = 8
						flag8 = true     // mark
						lastInfo = &info // leave this file process in the next batch
						// and we have to process files aleady in batch
					}
				} else {
					batch = append(batch, &info)
					if len(batch) < sBlock { // not filled
						continue
					}
				}

			}

			if len(batch) == 0 {
				break
			}

			b++

			prefix = fmt.Sprintf("[block #%03d]", b)

			wg0.Add(1)
			tokens0 <- 1
			go func(files []*UnikFileInfo, b int, prefix string) {
				var wg sync.WaitGroup
				tokens := make(chan int, opt.NumCPUs)

				var maxElements int64
				for _, info := range files {
					if maxElements < info.Kmers {
						maxElements = info.Kmers
					}
				}

				nBatchFiles := int((len(files) + 7) / 8)
				batchFiles := make([]string, 0, nBatchFiles)

				numSigs := CalcSignatureSize(uint64(maxElements), numHashes, fpr)
				var eFileSize float64
				eFileSize = 24
				for _, info := range files {
					eFileSize += float64(len(info.Name) + 1 + 8)
				}
				eFileSize += float64(numSigs * uint64(nBatchFiles))

				if opt.Verbose {
					log.Infof("%s #files: %d, max #k-mers: %d, #signatures: %d, file size: %s", prefix, len(files), maxElements, numSigs, bytesize.ByteSize(eFileSize))
				}

				var bb, jj int
				for ii := 0; ii < len(files); ii += 8 {
					if dryRun {
						continue
					}
					jj = ii + 8
					if jj > len(files) {
						jj = len(files)
					}
					wg.Add(1)
					tokens <- 1
					bb++

					outFile := filepath.Join(outDir, fmt.Sprintf("block%03d_batch%03d%s", b, bb, extIndex))
					batchFiles = append(batchFiles, outFile)

					go func(_files []*UnikFileInfo, bb int, maxElements int64, numSigs uint64, outFile string) {
						defer func() {
							wg.Done()
							<-tokens
						}()

						var outGzip bool
						if len(_files) > 1 {
							outGzip = true
						}
						outfh, gw, w, err := outStream(outFile, outGzip, opt.CompressionLevel)
						checkError(err)
						defer func() {
							outfh.Flush()
							if gw != nil {
								gw.Close()
							}
							w.Close()
						}()

						names := make([]string, 0, 8)
						sizes := make([]uint64, 0, 8)
						for _, info := range _files {
							names = append(names, info.Name)
							sizes = append(sizes, uint64(info.Kmers))
						}

						writer, err := index.NewWriter(outfh, k, canonical, uint8(numHashes), numSigs, names, sizes)
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
							var loc int

							tokensOpenFiles <- 1
							infh, r, _, err = inStream(info.Path)
							checkError(errors.Wrap(err, info.Path))

							reader, err = unikmer.NewReader(infh)
							checkError(errors.Wrap(err, info.Path))

							for {
								code, _, err = reader.ReadCodeWithTaxid()
								if err != nil {
									if err == io.EOF {
										break
									}
									checkError(errors.Wrap(err, info.Path))
								}

								for _, loc = range hashLocations(code, numHashes, numSigs) {
									sigs[loc] |= 1 << (7 - _k)
								}
							}

							r.Close()

							<-tokensOpenFiles
						}

						checkError(writer.WriteBatch(sigs, len(sigs)))
						if opt.Verbose {
							log.Infof("%s batch #%03d: %d signatures saved", prefix, bb, len(sigs))
						}

					}(files[ii:jj], bb, maxElements, numSigs, outFile)
				}

				wg.Wait()

				blockFile := filepath.Join(outDir, fmt.Sprintf("block%03d%s", b, extIndex))

				if !dryRun {
					if len(files) == 1 { // do not have to merge
						checkError(os.Rename(batchFiles[0], blockFile))
					} else {
						if opt.Verbose {
							log.Infof("%s merging %d index files", prefix, len(batchFiles))
						}

						checkError(MergeUnikIndex(opt, prefix, batchFiles, blockFile))

						if opt.Verbose {
							log.Infof("%s finished merging", prefix)
						}
					}
				}

				wg0.Done()
				ch <- filepath.Base(blockFile)
				chFileSize <- eFileSize
				<-tokens0
			}(batch, b, prefix)

			batch = make([]*UnikFileInfo, 0, sBlock)
		}

		wg0.Wait()
		close(ch)
		close(chFileSize)
		<-done
		<-doneFileSize

		if !dryRun {
			sort.Strings(indexFiles)
			dbInfo := NewUnikIndexDBInfo(int(index.Version), indexFiles)
			dbInfo.K = k
			dbInfo.Kmers = int(n)
			dbInfo.FPR = fpr
			dbInfo.BlockSize = sBlock0
			dbInfo.NumNames = len(names0)
			dbInfo.Names = names0
			dbInfo.Sizes = sizes0
			dbInfo.NumHashes = numHashes
			dbInfo.Canonical = canonical
			checkError(dbInfo.WriteTo(filepath.Join(outDir, dbInfoFile)))
		}

		// ------------------------------------------------------------------------------------

		if opt.Verbose {
			log.Infof("unikmer index database with %d k-mers saved to %s", n, outDir)
			log.Infof("#index files: %d, total file size: %s", len(indexFiles), bytesize.ByteSize(fileSize))
		}
	},
}

func init() {
	dbCmd.AddCommand(indexCmd)

	// indexCmd.Flags().StringP("out-prefix", "o", "-", `index file prefix ("-" for stdout)`)
	indexCmd.Flags().StringP("out-dir", "O", "unikmer-db", `output directory`)
	indexCmd.Flags().Float64P("false-positive-rate", "f", 0.3, `false positive rate of single bloom filter`)
	indexCmd.Flags().IntP("num-hash", "n", 1, `number of hashes`)
	indexCmd.Flags().IntP("block-size", "b", 0, `block size, better be multiple of 64. default: sqrt(#.files)`)
	indexCmd.Flags().StringP("block-max-kmers-t1", "m", "20M", `if kmers of single .unik file exceeds this threshold, we creat change block size to 8. unit supported: K, M, G`)
	indexCmd.Flags().StringP("block-max-kmers-t2", "M", "200M", `if kmers of single .unik file exceeds this threshold, we creat individual index for this file. unit supported: K, M, G`)

	indexCmd.Flags().BoolP("force", "", false, "overwrite tmp dir")
	indexCmd.Flags().StringP("name-regexp", "r", "", "regular expression for extract name from .unik file name. if not given, base name are saved")
	indexCmd.Flags().IntP("max-open-files", "F", 256, "maximum number of opened files")
	indexCmd.Flags().BoolP("dry-run", "", false, "dry run, useful to adjust parameters")
}
