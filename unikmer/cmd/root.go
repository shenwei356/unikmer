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
	"fmt"
	"os"
	"runtime"

	"github.com/klauspost/compress/flate"
	homedir "github.com/mitchellh/go-homedir"
	"github.com/spf13/cobra"
)

// RootCmd represents the base command when called without any subcommands
var RootCmd = &cobra.Command{
	Use:   "unikmer",
	Short: "Unique-Kmer Toolkit",
	Long: fmt.Sprintf(`unikmer - Unique-Kmer Toolkit

A command-line toolkit providing functions including counting, format
convertion, set operations and searching of small k-mers (k <= 32)
optional with Taxids but without frequency information.

K-mers (k <= 32) are encoded into 'uint64', stored in builtin 'map' of
golang in RAM, and serialized in binary file with extension '.unik'.

Taxids can be assigned when counting k-mers from genome sequences,
and LCA (Lowest Common Ancestor) will be computed during opertions
including computing union, intersecton, set difference, unique and
repeated k-mers.

Version: v%s

Author: Wei Shen <shenwei356@gmail.com>

Documents  : https://shenwei356.github.io/unikmer
Source code: https://github.com/shenwei356/unikmer

Dataset (optional):

  Some commands need taxonomy nodes file from e.g., NCBI Taxonomy database,
  you can extract "nodes.dmp" from link below into ~/.unikmer/ ,
  ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz , 
  or some other directory, and later you can refer to using flag --data-dir,
  or environment variable UNIKMER_DB.

  For GTDB, use https://github.com/nick-youngblut/gtdb_to_taxdump 
  for taxonomy convertion.

  Note that Taxids are represented using uint32 and stored in 4 or less bytes,
  all taxids should be in range of [1, %d]

`, VERSION, maxUint32),
}

// Execute adds all child commands to the root command sets flags appropriately.
// This is called by main.main(). It only needs to happen once to the rootCmd.
func Execute() {
	if err := RootCmd.Execute(); err != nil {
		fmt.Println(err)
		os.Exit(-1)
	}
}

var defaulDataDir string

func init() {
	var err error
	defaulDataDir, err = homedir.Expand("~/.unikmer/")
	checkError(err)

	defaultThreads := runtime.NumCPU()
	if defaultThreads > 2 {
		defaultThreads = 2
	}

	RootCmd.PersistentFlags().IntP("threads", "j", defaultThreads, "number of CPUs to use. (default value: 1 for single-CPU PC, 2 for others)")
	RootCmd.PersistentFlags().BoolP("verbose", "", false, "print verbose information")
	RootCmd.PersistentFlags().BoolP("no-compress", "C", false, "do not compress binary file (not recommended)")
	RootCmd.PersistentFlags().IntP("compression-level", "L", flate.DefaultCompression, "compression level")
	RootCmd.PersistentFlags().BoolP("compact", "c", false, "write more compact binary file with little loss of speed")
	RootCmd.PersistentFlags().StringP("infile-list", "i", "", "file of input files list (one file per line), if given, they are appended to files from cli arguments")

	RootCmd.PersistentFlags().Uint32P("max-taxid", "", 1<<32-1, "for smaller taxids, we can use less space to store taxids. default value is 1<<32-1, that's enough for NCBI Taxonomy taxids")
	RootCmd.PersistentFlags().BoolP("ignore-taxid", "I", false, "ignore taxonomy information")
	RootCmd.PersistentFlags().StringP("data-dir", "", defaulDataDir, "directory containing nodes.dmp and names.dmp")
	// RootCmd.PersistentFlags().BoolP("cache-lca", "", false, "cache LCA queries")
}

const helpSort = "sort k-mers, this significantly reduce file size. You can even disable gzip compression by flag -C/--no-compress. This flag overwrites global flag -c/--compact"
