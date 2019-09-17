// Copyright © 2018-2019 Wei Shen <shenwei356@gmail.com>
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

	"github.com/shenwei356/util/pathutil"

	"github.com/shenwei356/unikmer"
	"github.com/spf13/cobra"
)

// mergeCmd represents
var mergeCmd = &cobra.Command{
	Use:   "merge",
	Short: "merge from sorted chunk files",
	Long: `merge from sorted chunk files

`,
	Run: func(cmd *cobra.Command, args []string) {
		opt := getOptions(cmd)
		runtime.GOMAXPROCS(opt.NumCPUs)

		outFile0 := getFlagString(cmd, "out-prefix")
		unique := getFlagBool(cmd, "unique")
		repeated := getFlagBool(cmd, "repeated")
		maxOpenFiles := getFlagPositiveInt(cmd, "max-open-files")

		var err error

		var files []string
		infileList := getFlagString(cmd, "infile-list")
		if infileList != "" {
			files, err = getListFromFile(infileList)
			checkError(err)
		} else {
			files = getFileList(args)
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
			log.Infof("checking chunk files")
		}
		var infh *bufio.Reader
		var r *os.File
		var reader *unikmer.Reader
		var nUnequalK, nNotConsC, nNotSorted int
		var k int = -1
		var canonical bool
		var mode uint32

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

				reader, err = unikmer.NewReader(infh)
				checkError(err)

				if k == -1 { // first file
					k = reader.K
					canonical = reader.Flag&unikmer.UNIK_CANONICAL > 0

					if opt.Compact {
						mode |= unikmer.UNIK_COMPACT
					}
					if canonical {
						mode |= unikmer.UNIK_CANONICAL
					}
					mode |= unikmer.UNIK_SORTED
				} else if k != reader.K {
					log.Errorf("K (%d) of binary file '%s' not equal to previous K (%d)", reader.K, file, k)
					nUnequalK++
				} else if (reader.Flag&unikmer.UNIK_CANONICAL > 0) != canonical {
					log.Errorf(`'canonical' flags not consistent, please check with "unikmer stats"`)
					nNotConsC++
				} else if reader.Flag&unikmer.UNIK_SORTED == 0 { // not sorted
					log.Errorf("chunk file not sorted: %s", file)
					nNotSorted++
				}
			}()
		}
		if nNotSorted > 0 || nUnequalK > 0 || nNotConsC > 0 {
			checkError(fmt.Errorf("please check chunk files: %d with different K, %d with different canonical flag, %d not sorted", nUnequalK, nNotConsC, nNotSorted))
		}

		files = _files

		if len(files) == 0 {
			log.Warningf("no valid chunk files given")
			return
		}
		if opt.Verbose {
			log.Infof("merging from %d chunk files", len(files))
		}

		// merge

		outFile := outFile0
		if !isStdout(outFile) {
			outFile += extDataFile
		}

		if len(files) < maxOpenFiles {
			n, _ := mergeChunksFile(opt, files, outFile, k, mode, unique, repeated)

			if opt.Verbose {
				log.Infof("%d k-mers saved to %s", n, outFile)
			}
			return
		}

		if opt.Verbose {
			log.Infof("two-pass merge performing")
		}

		if maxOpenFiles > len(files)*len(files) {
			log.Warningf("are you sure for merging from %d files?", len(files)*len(files))
			log.Warningf("if the files are of small size, you may use 'unikmer sort -m' instead")
		}

		tmpDir := getFlagString(cmd, "tmp-dir")
		if isStdout(outFile0) {
			tmpDir = filepath.Join(tmpDir, "stdout.tmp")
		} else {
			tmpDir = filepath.Join(tmpDir, filepath.Base(outFile0)+".tmp")
		}
		err = os.MkdirAll(tmpDir, 0755)
		if err != nil {
			checkError(fmt.Errorf("fail to create temp directory: %s", tmpDir))
		}

		tmpFiles := make([]string, 0, 10)
		iTmpFile := 0
		_files = make([]string, 0, maxOpenFiles)
		for _, file := range files {
			_files = append(_files, file)
			if len(_files) == maxOpenFiles {
				iTmpFile++
				outFile1 := chunkFileName(tmpDir, iTmpFile)

				n, _ := mergeChunksFile(opt, _files, outFile1, k, mode, false, false)
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

			n, _ := mergeChunksFile(opt, _files, outFile1, k, mode, false, false)
			if opt.Verbose {
				log.Infof("%d k-mers saved to tmp file: %s", n, outFile1)
			}
			tmpFiles = append(tmpFiles, outFile1)
		}

		n, _ := mergeChunksFile(opt, tmpFiles, outFile, k, mode, unique, repeated)

		if opt.Verbose {
			log.Infof("%d k-mers saved to %s", n, outFile)
		}

		// cleanning

		if opt.Verbose {
			log.Infof("removing %d intermediate files", len(tmpFiles))
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
			checkError(fmt.Errorf("fail to remove temp directory: %s", tmpDir))
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

	mergeCmd.Flags().IntP("max-open-files", "M", 300, `max number of open files`)
	mergeCmd.Flags().StringP("tmp-dir", "t", "./", `directory for intermediate files`)
}
