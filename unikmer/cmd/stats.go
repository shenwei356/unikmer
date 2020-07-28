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
	"sort"
	"strconv"
	"strings"
	"sync"

	"github.com/cznic/sortutil"
	humanize "github.com/dustin/go-humanize"
	"github.com/shenwei356/unikmer"
	"github.com/spf13/cobra"
	prettytable "github.com/tatsushid/go-prettytable"
)

// statCmd represents
var statCmd = &cobra.Command{
	Use:   "stats",
	Short: "Statistics of binary files",
	Long: `Statistics of binary files

Tips:
  1. For lots of small files (especially on SDD), use big value of '-j' to
     parallelize counting.

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
		all := getFlagBool(cmd, "all")
		tabular := getFlagBool(cmd, "tabular")
		skipErr := getFlagBool(cmd, "skip-err")
		sTrue := getFlagString(cmd, "symbol-true")
		sFalse := getFlagString(cmd, "symbol-false")
		basename := getFlagBool(cmd, "basename")

		if sTrue == sFalse {
			checkError(fmt.Errorf("values of -/--symbol-true and -F/--symbol--false should be different"))
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

		// tabular output
		if tabular {
			colnames := []string{
				"file",
				"k",
				"gzipped",
				"compact",
				"canonical",
				"sorted",
				"include-taxid",
				"global-taxid",
			}
			if all {
				colnames = append(colnames, []string{"number"}...)
			}
			outfh.WriteString(strings.Join(colnames, "\t") + "\n")
			outfh.Flush()
		}

		ch := make(chan statInfo, opt.NumCPUs)
		statInfos := make([]statInfo, 0, 256)

		cancel := make(chan struct{})

		done := make(chan int)

		go func() {
			var id uint64 = 1 // for keepping order
			buf := make(map[uint64]statInfo)

			for info := range ch {
				if info.err != nil {
					if skipErr {
						log.Warningf("%s: %s", info.file, info.err)
						continue
					} else {
						log.Errorf("%s: %s", info.file, info.err)
						close(cancel)
						break
					}
				}

				if id == info.id { // right the one
					if !tabular {
						statInfos = append(statInfos, info)
					} else {
						if !all {
							outfh.WriteString(fmt.Sprintf("%s\t%v\t%v\t%v\t%v\t%v\t%v\t%s\n",
								info.file,
								info.k,
								boolStr(sTrue, sFalse, info.gzipped),
								boolStr(sTrue, sFalse, info.compact),
								boolStr(sTrue, sFalse, info.canonical),
								boolStr(sTrue, sFalse, info.sorted),
								boolStr(sTrue, sFalse, info.includeTaxid),
								info.globalTaxid,
							))
						} else {
							outfh.WriteString(fmt.Sprintf("%s\t%v\t%v\t%v\t%v\t%v\t%v\t%s\t%d\n",
								info.file,
								info.k,
								boolStr(sTrue, sFalse, info.gzipped),
								boolStr(sTrue, sFalse, info.compact),
								boolStr(sTrue, sFalse, info.canonical),
								boolStr(sTrue, sFalse, info.sorted),
								boolStr(sTrue, sFalse, info.includeTaxid),
								info.globalTaxid,
								info.number,
							))
						}
						outfh.Flush()
					}
					id++
				} else { // check bufferd result
					for true {
						if info1, ok := buf[id]; ok {
							if !tabular {
								statInfos = append(statInfos, info1)
							} else {
								if !all {
									outfh.WriteString(fmt.Sprintf("%s\t%v\t%v\t%v\t%v\t%v\t%v\t%s\n",
										info.file,
										info.k,
										boolStr(sTrue, sFalse, info.gzipped),
										boolStr(sTrue, sFalse, info.compact),
										boolStr(sTrue, sFalse, info.canonical),
										boolStr(sTrue, sFalse, info.sorted),
										boolStr(sTrue, sFalse, info.includeTaxid),
										info.globalTaxid,
									))
								} else {
									outfh.WriteString(fmt.Sprintf("%s\t%v\t%v\t%v\t%v\t%v\t%v\t%s\t%d\n",
										info.file,
										info.k,
										boolStr(sTrue, sFalse, info.gzipped),
										boolStr(sTrue, sFalse, info.compact),
										boolStr(sTrue, sFalse, info.canonical),
										boolStr(sTrue, sFalse, info.sorted),
										boolStr(sTrue, sFalse, info.includeTaxid),
										info.globalTaxid,
										info.number,
									))
								}
								outfh.Flush()
							}

							delete(buf, info1.id)
							id++
						} else {
							break
						}
					}
					buf[info.id] = info
				}
			}

			if len(buf) > 0 {
				ids := make(sortutil.Uint64Slice, len(buf))
				i := 0
				for id := range buf {
					ids[i] = id
					i++
				}
				sort.Sort(ids)
				for _, id := range ids {
					info := buf[id]
					if !tabular {
						statInfos = append(statInfos, info)
					} else {
						if !all {
							outfh.WriteString(fmt.Sprintf("%s\t%v\t%v\t%v\t%v\t%v\t%v\t%s\n",
								info.file,
								info.k,
								boolStr(sTrue, sFalse, info.gzipped),
								boolStr(sTrue, sFalse, info.compact),
								boolStr(sTrue, sFalse, info.canonical),
								boolStr(sTrue, sFalse, info.sorted),
								boolStr(sTrue, sFalse, info.includeTaxid),
								info.globalTaxid,
							))
						} else {
							outfh.WriteString(fmt.Sprintf("%s\t%v\t%v\t%v\t%v\t%v\t%v\t%s\t%d\n",
								info.file,
								info.k,
								boolStr(sTrue, sFalse, info.gzipped),
								boolStr(sTrue, sFalse, info.compact),
								boolStr(sTrue, sFalse, info.canonical),
								boolStr(sTrue, sFalse, info.sorted),
								boolStr(sTrue, sFalse, info.includeTaxid),
								info.globalTaxid,
								info.number,
							))
						}
						outfh.Flush()
					}
				}
			}

			done <- 1
		}()

		chFile := make(chan string, opt.NumCPUs)
		doneSendFile := make(chan int)
		go func() {
			for _, file := range files {
				select {
				case <-cancel:
					break
				default:
				}
				chFile <- file
			}
			close(chFile)
			doneSendFile <- 1
		}()

		var wg sync.WaitGroup
		token := make(chan int, opt.NumCPUs)
		var id uint64

		for file := range chFile {
			select {
			case <-cancel:
				break
			default:
			}

			token <- 1
			wg.Add(1)
			id++
			go func(file string, id uint64) {
				defer func() {
					wg.Done()
					<-token
				}()

				var infh *bufio.Reader
				var r *os.File
				var reader *unikmer.Reader
				var gzipped bool
				var n int64
				var globalTaxid string

				infh, r, gzipped, err = inStream(file)
				if err != nil {
					select {
					case <-cancel:
						return
					default:
					}
					if basename {
						file = filepath.Base(file)
					}
					ch <- statInfo{file: file, err: err, id: id}
					return
				}
				defer r.Close()

				reader, err = unikmer.NewReader(infh)
				checkError(err)
				if err != nil {
					select {
					case <-cancel:
						return
					default:
					}

					if basename {
						file = filepath.Base(file)
					}
					ch <- statInfo{file: file, err: err, id: id}
					return
				}

				n = 0
				if all {
					if reader.IsSorted() && reader.Number >= 0 {
						n = reader.Number
					} else {
						for {
							_, _, err = reader.ReadCodeWithTaxid()
							if err != nil {
								if err == io.EOF {
									break
								}
								checkError(err)
							}

							n++
						}
					}
				}
				if basename {
					file = filepath.Base(file)
				}
				if reader.GetGlobalTaxid() > 0 {
					globalTaxid = strconv.FormatUint(uint64(reader.GetGlobalTaxid()), 10)
				} else {
					globalTaxid = ""
				}
				ch <- statInfo{
					file:         file,
					k:            reader.K,
					gzipped:      gzipped,
					compact:      reader.IsCompact(),
					canonical:    reader.IsCanonical(),
					sorted:       reader.IsSorted(),
					includeTaxid: reader.IsIncludeTaxid(),
					globalTaxid:  globalTaxid,
					number:       n,

					err: nil,
					id:  id,
				}

			}(file, id)
		}

		<-doneSendFile
		wg.Wait()
		close(ch)
		<-done

		select {
		case <-cancel:
			return
		default:
		}

		if tabular {
			return
		}

		// format output
		columns := []prettytable.Column{
			{Header: "file"},
			{Header: "k", AlignRight: true},
			{Header: "gzipped"},
			{Header: "compact"},
			{Header: "canonical"},
			{Header: "sorted"},
			{Header: "include-taxid"},
			{Header: "global-taxid"},
		}
		if all {
			columns = append(columns, []prettytable.Column{
				{Header: "number", AlignRight: true},
			}...)
		}
		tbl, err := prettytable.NewTable(columns...)

		checkError(err)
		tbl.Separator = "  "

		for _, info := range statInfos {
			if !all {
				tbl.AddRow(
					info.file,
					info.k,
					boolStr(sTrue, sFalse, info.gzipped),
					boolStr(sTrue, sFalse, info.compact),
					boolStr(sTrue, sFalse, info.canonical),
					boolStr(sTrue, sFalse, info.sorted),
					boolStr(sTrue, sFalse, info.includeTaxid),
					info.globalTaxid,
				)
			} else {
				tbl.AddRow(
					info.file,
					info.k,
					boolStr(sTrue, sFalse, info.gzipped),
					boolStr(sTrue, sFalse, info.compact),
					boolStr(sTrue, sFalse, info.canonical),
					boolStr(sTrue, sFalse, info.sorted),
					boolStr(sTrue, sFalse, info.includeTaxid),
					info.globalTaxid,
					humanize.Comma(info.number),
				)
			}
		}
		outfh.Write(tbl.Bytes())
	},
}

type statInfo struct {
	file         string
	k            int
	gzipped      bool
	compact      bool
	canonical    bool
	sorted       bool
	includeTaxid bool
	globalTaxid  string
	number       int64

	err error
	id  uint64
}

func init() {
	RootCmd.AddCommand(statCmd)

	statCmd.Flags().StringP("out-file", "o", "-", `out file ("-" for stdout, suffix .gz for gzipped out)`)
	statCmd.Flags().BoolP("all", "a", false, "all information, including number of k-mers")
	statCmd.Flags().BoolP("tabular", "T", false, "output in machine-friendly tabular format")
	statCmd.Flags().BoolP("skip-err", "e", false, "skip error, only show warning message")
	statCmd.Flags().StringP("symbol-true", "", "✓", "smybol for true")
	statCmd.Flags().StringP("symbol-false", "", "✕", "smybol for false")
	statCmd.Flags().BoolP("basename", "b", false, "only output basename of files")
}

func boolStr(sTrue, sFalse string, v bool) string {
	if v {
		return sTrue
	}
	return sFalse
}
