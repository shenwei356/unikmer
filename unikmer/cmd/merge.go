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
	"io/ioutil"
	"os"
	"path/filepath"
	"regexp"
	"runtime"

	"github.com/pkg/errors"
	"github.com/shenwei356/util/pathutil"

	"github.com/shenwei356/unikmer"
	"github.com/spf13/cobra"
)

var mergeCmd = &cobra.Command{
	Use:   "merge",
	Short: "Merge k-mers from sorted chunk files",
	Long: `Merge k-mers from sorted chunk files

Attentions:
  1. The 'canonical' flags of all files should be consistent.
  2. Input files should ALL have or don't have taxid information.
  3. Input files should be sorted.
  
Tips:
  1. If you don't need to compute unique or repeated k-mers, 
     use 'unikmer concat -s', which is faster.
`,
	Run: func(cmd *cobra.Command, args []string) {
		opt := getOptions(cmd)
		runtime.GOMAXPROCS(opt.NumCPUs)

		outFile0 := getFlagString(cmd, "out-prefix")
		unique := getFlagBool(cmd, "unique")
		repeated := getFlagBool(cmd, "repeated")
		maxOpenFiles := getFlagPositiveInt(cmd, "max-open-files")
		keepTmpDir := getFlagBool(cmd, "keep-tmp-dir")
		force := getFlagBool(cmd, "force")

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

		// read files from directory
		if getFlagBool(cmd, "is-dir") {
			paternChunkFile := getFlagString(cmd, "pattern")
			reChunFile, err := regexp.Compile(paternChunkFile)
			if err != nil {
				checkError(fmt.Errorf("fail to compile pattern of chunk file: %s", paternChunkFile))
			}

			_files := make([]string, 0, 10)
			var existed bool
			var list []os.FileInfo
			var filename string
			var N, n int
			for _, dir := range files {
				if isStdin(dir) {
					dir = "./"
				}
				N++
				n = 0
				existed, err = pathutil.DirExists(dir)
				if err != nil {
					checkError(fmt.Errorf("check given dir '%s': %s", dir, err))
				}
				if !existed {
					log.Warning("skip unexisted dir: %s", dir)
				}

				list, err = ioutil.ReadDir(dir)
				if err != nil {
					checkError(fmt.Errorf("check given directory '%s': %s", dir, err))
				}

				for _, file := range list {
					filename = file.Name()

					if filename[0] == '.' || file.IsDir() || !reChunFile.MatchString(filename) {
						continue
					}

					_files = append(_files, filepath.Join(dir, filename))
					n++
				}
				if opt.Verbose {
					log.Infof("%d chunk files found in dir: %s", n, dir)
				}
			}
			files = _files

			if n == 0 {
				log.Warningf("%d chunk files found in %d dir(s)", n, N)
				return
			} else if opt.Verbose {
				log.Infof("%d chunk files found in %d dir(s)", n, N)
			}

		}

		// checking files
		if opt.Verbose {
			log.Info()
			log.Infof("======= Stage 1: checking chunk files =======")
		}
		var infh *bufio.Reader
		var r *os.File
		var reader0 *unikmer.Reader
		var k int = -1
		var canonical bool
		var hashed bool
		var hasTaxid bool
		var mode uint32
		var taxondb *unikmer.Taxonomy

		_files := make([]string, 0, len(files))
		for _, file := range files {
			if isStdout(file) {
				log.Warningf("skip stdin")
				continue
			}
			_files = append(_files, file)
			func() {
				infh, r, _, err = inStream(file)
				checkError(err)
				defer r.Close()

				reader, err := unikmer.NewReader(infh)
				checkError(errors.Wrap(err, file))

				if !reader.IsSorted() {
					checkError(fmt.Errorf("input files should be sorted"))
				}

				if k == -1 { // first file
					reader0 = reader
					k = reader.K
					canonical = reader.IsCanonical()
					hashed = reader.IsHashed()
					hasTaxid = !opt.IgnoreTaxid && reader.HasTaxidInfo()

					if canonical {
						mode |= unikmer.UNIK_CANONICAL
					}
					if hasTaxid {
						mode |= unikmer.UNIK_INCLUDETAXID
					}
					if hashed {
						mode |= unikmer.UNIK_HASHED
					}
					mode |= unikmer.UNIK_SORTED

					if hasTaxid {
						if opt.Verbose {
							log.Infof("taxids found in file: %s", file)
						}
						taxondb = loadTaxonomy(opt, false)
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

			}()
		}

		files = _files

		if len(files) == 0 {
			log.Warningf("no valid chunk files given")
			return
		}
		if opt.Verbose {
			log.Infof("checking passed")
		}

		// merge

		outFile := outFile0
		if !isStdout(outFile) {
			outFile += extDataFile
		}

		if len(files) < maxOpenFiles {
			if opt.Verbose {
				log.Info()
				log.Infof("======= Stage 2: merging from %d chunks =======", len(files))
			}
			n, _ := mergeChunksFile(opt, taxondb, files, outFile, k, mode, unique, repeated, true)

			if opt.Verbose {
				log.Infof("%d k-mers saved to %s", n, outFile)
			}
			return
		}

		if opt.Verbose {
			log.Info()
			log.Infof("======= Stage 2: merging from %d chunks (round: 1/2) =======", len(files))
		}

		// if maxOpenFiles > len(files)*len(files) {
		// 	log.Warningf("are you sure for merging from %d files?", len(files)*len(files))
		// 	log.Warningf("if the files are of small size, you may use 'unikmer sort -m' instead")
		// }

		tmpDir := getFlagString(cmd, "tmp-dir")
		if isStdout(outFile0) {
			tmpDir = filepath.Join(tmpDir, "stdout.tmp")
		} else {
			tmpDir = filepath.Join(tmpDir, filepath.Base(outFile0)+".tmp")
		}

		existed, err := pathutil.DirExists(tmpDir)
		checkError(errors.Wrap(err, tmpDir))
		if existed {
			empty, err := pathutil.IsEmpty(tmpDir)
			checkError(errors.Wrap(err, tmpDir))
			if !empty {
				if force {
					checkError(os.RemoveAll(tmpDir))
				} else {
					checkError(fmt.Errorf("tmp dir not empty: %s, choose another one or use --force to overwrite", tmpDir))
				}
			} else {
				checkError(os.RemoveAll(tmpDir))
			}
		}
		checkError(os.MkdirAll(tmpDir, 0777))

		tmpFiles := make([]string, 0, 10)
		iTmpFile := 0
		_files = make([]string, 0, maxOpenFiles)
		for _, file := range files {
			_files = append(_files, file)
			if len(_files) == maxOpenFiles {
				iTmpFile++
				outFile1 := chunkFileName(tmpDir, iTmpFile)

				if opt.Verbose {
					log.Infof("[chunk %d] merging k-mers from %d tmp files", iTmpFile, len(_files))
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
				log.Infof("[chunk %d] merging k-mers from %d tmp files", iTmpFile, len(_files))
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
		n, _ := mergeChunksFile(opt, taxondb, tmpFiles, outFile, k, mode, unique, repeated, true)

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
		for _, file := range tmpFiles {
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

	},
}

func init() {
	RootCmd.AddCommand(mergeCmd)

	mergeCmd.Flags().BoolP("is-dir", "D", false, "intput files are directory containing chunk files")
	mergeCmd.Flags().StringP("pattern", "p", `^chunk_\d+\.unik$`, `chunk file pattern (regular expression)`)

	mergeCmd.Flags().StringP("out-prefix", "o", "-", `out file prefix ("-" for stdout)`)
	mergeCmd.Flags().BoolP("unique", "u", false, `remove duplicated k-mers`)
	mergeCmd.Flags().BoolP("repeated", "d", false, `only print duplicate k-mers`)

	mergeCmd.Flags().IntP("max-open-files", "M", 400, `max number of open files`)
	mergeCmd.Flags().StringP("tmp-dir", "t", "./", `directory for intermediate files`)
	mergeCmd.Flags().BoolP("keep-tmp-dir", "k", false, `keep tmp dir`)
	mergeCmd.Flags().BoolP("force", "", false, "overwrite tmp dir")
}
