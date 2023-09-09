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
	Short: "Toolkit for k-mer with taxonomic information",
	Long: fmt.Sprintf(`unikmer - Toolkit for k-mer with taxonomic information

unikmer is a toolkit for nucleic acid k-mer analysis, providing functions
including set operation on k-mers optional with TaxIds but without count
information.

K-mers are either encoded (k<=32) or hashed (k<=64) into 'uint64',
and serialized in binary file with the extension '.unik'.

TaxIds can be assigned when counting k-mers from genome sequences,
and LCA (Lowest Common Ancestor) is computed during set opertions
including computing union, intersection, set difference, unique and
repeated k-mers.

Version: v%s

Author: Wei Shen <shenwei356@gmail.com>

Documents  : https://bioinf.shenwei.me/unikmer
Source code: https://github.com/shenwei356/unikmer

Dataset (optional):

  Manipulating k-mers with TaxIds needs taxonomy file from e.g., 
  NCBI Taxonomy database, please extract "nodes.dmp", "names.dmp",
  "delnodes.dmp" and "merged.dmp" from link below into ~/.unikmer/ ,
  ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz , 
  or some other directory, and later you can refer to using flag
  --data-dir or environment variable UNIKMER_DB.

  For GTDB, use 'taxonkit create-taxdump' to create NCBI-style
  taxonomy dump files, or download from:
    https://github.com/shenwei356/gtdb-taxonomy

  Note that TaxIds are represented using uint32 and stored in 4 or
  less bytes, all TaxIds should be in the range of [1, %d].

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

var defaultDataDir string

func init() {
	var err error
	defaultDataDir, err = homedir.Expand("~/.unikmer/")
	checkError(err)

	defaultThreads := runtime.NumCPU()
	if defaultThreads > 4 {
		defaultThreads = 4
	}

	RootCmd.PersistentFlags().IntP("threads", "j", defaultThreads, "number of CPUs to use")
	RootCmd.PersistentFlags().BoolP("verbose", "", false, "print verbose information")
	RootCmd.PersistentFlags().BoolP("no-compress", "C", false, "do not compress binary file (not recommended)")
	RootCmd.PersistentFlags().IntP("compression-level", "", flate.DefaultCompression, "compression level")
	RootCmd.PersistentFlags().BoolP("compact", "c", false, "write compact binary file with little loss of speed")
	RootCmd.PersistentFlags().StringP("infile-list", "i", "", "file of input files list (one file per line), if given, they are appended to files from cli arguments")

	RootCmd.PersistentFlags().Uint32P("max-taxid", "", 1<<32-1, "for smaller TaxIds, we can use less space to store TaxIds. default value is 1<<32-1, that's enough for NCBI Taxonomy TaxIds")
	RootCmd.PersistentFlags().BoolP("ignore-taxid", "I", false, "ignore taxonomy information")
	RootCmd.PersistentFlags().StringP("data-dir", "", defaultDataDir, "directory containing NCBI Taxonomy files, including nodes.dmp, names.dmp, merged.dmp and delnodes.dmp")

	RootCmd.PersistentFlags().BoolP("nocheck-file", "", false, "do not check binary file, when using process substitution or named pipe")

	RootCmd.CompletionOptions.DisableDefaultCmd = true
	RootCmd.SetHelpCommand(&cobra.Command{Hidden: true})

	RootCmd.SetUsageTemplate(usageTemplate(""))
}

const helpSort = "sort k-mers, this significantly reduce file size for k<=25. This flag overides global flag -c/--compact"

func usageTemplate(s string) string {
	return fmt.Sprintf(`Usage:{{if .Runnable}}
  {{.UseLine}}{{end}}{{if .HasAvailableSubCommands}}
  {{.CommandPath}} [command]{{end}} %s{{if gt (len .Aliases) 0}}

Aliases:
  {{.NameAndAliases}}{{end}}{{if .HasExample}}

Examples:
{{.Example}}{{end}}{{if .HasAvailableSubCommands}}

Available Commands:{{range .Commands}}{{if (or .IsAvailableCommand (eq .Name "help"))}}
  {{rpad .Name .NamePadding }} {{.Short}}{{end}}{{end}}{{end}}{{if .HasAvailableLocalFlags}}

Flags:
{{.LocalFlags.FlagUsagesWrapped 110 | trimTrailingWhitespaces}}{{end}}{{if .HasAvailableInheritedFlags}}

Global Flags:
{{.InheritedFlags.FlagUsagesWrapped 110 | trimTrailingWhitespaces}}{{end}}{{if .HasHelpSubCommands}}

Additional help topics:{{range .Commands}}{{if .IsAdditionalHelpTopicCommand}}
  {{rpad .CommandPath .CommandPathPadding}} {{.Short}}{{end}}{{end}}{{end}}{{if .HasAvailableSubCommands}}

Use "{{.CommandPath}} [command] --help" for more information about a command.{{end}}
`, s)
}
