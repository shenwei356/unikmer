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
	"fmt"
	"os"
	"runtime"

	"github.com/spf13/cobra"
)

// RootCmd represents the base command when called without any subcommands
var RootCmd = &cobra.Command{
	Use:   "unikmer",
	Short: "Unique-Kmer Toolkit",
	Long: fmt.Sprintf(`unikmer - Unique-Kmer Toolkit

A command-line toolkit providing functions including counting, format
convertion, set operations and searching on unique Kmers (k <= 32) while
NOT recording Kmer frequencies.

Kmers (k <= 32) are encoded into 'uint64', stored in builtin 'map' of golang
in RAM, and serialized in binary format.

Version: %s

Author: Wei Shen <shenwei356@gmail.com>

Documents  : https://shenwei356.github.io/unikmer
Source code: https://github.com/shenwei356/unikmer


`, VERSION),
}

// Execute adds all child commands to the root command sets flags appropriately.
// This is called by main.main(). It only needs to happen once to the rootCmd.
func Execute() {
	if err := RootCmd.Execute(); err != nil {
		fmt.Println(err)
		os.Exit(-1)
	}
}

func init() {
	defaultThreads := runtime.NumCPU()
	if defaultThreads > 2 {
		defaultThreads = 2
	}

	RootCmd.PersistentFlags().IntP("threads", "j", defaultThreads, "number of CPUs to use. (default value: 1 for single-CPU PC, 2 for others)")
	RootCmd.PersistentFlags().BoolP("verbose", "", false, "print verbose information")
	RootCmd.PersistentFlags().BoolP("no-compress", "C", false, "do not compress binary file (not recommended)")
	RootCmd.PersistentFlags().BoolP("compact", "c", false, "write more compact binary file with little loss of speed")
	RootCmd.PersistentFlags().StringP("infile-list", "i", "", "file of input files list (one file per line), if given, files from cli arguments are ignored")
}

const helpSort = "sort Kmers, this significantly reduce file size. You can even disable gzip compression by flag -C/--no-compress. This flag overwrites global option -c/--compact"
