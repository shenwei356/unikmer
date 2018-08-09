# unikmer

unikmer (unique Kmer) is a golang package and a command-line toolkit for
manipulating [Kmers](https://en.wikipedia.org/wiki/K-mer) while NOT recording
Kmer frequencies.

## The package

[![GoDoc](https://godoc.org/github.com/shenwei356/unikmer?status.svg)](https://godoc.org/github.com/shenwei356/unikmer)
[![Go Report Card](https://goreportcard.com/badge/github.com/shenwei356/unikmer)](https://goreportcard.com/report/github.com/shenwei356/unikmer)

The unikmer package provides basic manipulations of unique Kmers (NOT including
Kmer frequencies) and its binary file.

### Installation

    go get -u github.com/shenwei356/unikmer

### Benchmark

    $ go test . -bench=Bench*
    goos: linux
    goarch: amd64
    pkg: github.com/shenwei356/unikmer
    BenchmarkEncodeK32-4            20000000                98.1 ns/op
    BenchmarkDecodeK32-4            20000000                 102 ns/op
    BenchmarkRevK32-4               20000000                64.2 ns/op
    BenchmarkCompK32-4              20000000                54.8 ns/op
    BenchmarkRevCompK32-4           10000000                 116 ns/op


## The toolkit

`unikmer` is a command-line toolkit provides some functions including counting,
format convertion, set operations and searching on unique Kmers.

### Installation

1. Download [binary files](https://github.com/shenwei356/unikmer/releases).

1. Bioconda (not available now)

        conda install unikmer

### Commands

1. Counting

        count           count Kmer from FASTA/Q sequences
        subset          extract smaller Kmers from binary file

1. Format convertion

        view            read and output binary format to plain text
        dump            convert plain Kmer text to binary format

1. Set operations

        concat          concatenate multiple binary files
        diff            set difference of multiple binary files
        inter           intersection of multiple binary files
        union           union of multiple binary files

1. Searching

        grep            search Kmer from binary file

1. Misc

        genautocomplete generate shell autocompletion script
        help            Help about any command
        version         print version information and check for update

### Quick Start

    # counting
    $ time unikmer count -k 31 Ecoli-MG1655.fasta.gz -o Ecoli-MG1655.fasta.gz
    real    0m5.209s
    user    0m6.864s
    sys     0m0.169s

    $ ls -lh Ecoli-MG1655.fasta.gz*
    -rw-rw-r--. 1 shenwei shenwei 1.4M Aug  9 23:19 Ecoli-MG1655.fasta.gz
    -rw-rw-r--. 1 shenwei shenwei  23M Aug  9 23:29 Ecoli-MG1655.fasta.gz.unik


    # view
    $ unikmer view Ecoli-MG1655.fasta.gz.unik | head -n 3
    AGCTTTTCATTCTGACTGCAACGGGCAATAT
    GCTTTTCATTCTGACTGCAACGGGCAATATG
    CTTTTCATTCTGACTGCAACGGGCAATATGT

    $ unikmer view Ecoli-MG1655.fasta.gz.unik | wc -l
    9108538


    # union
    $ unikmer union Ecoli-MG1655.fasta.gz.unik Ecoli-IAI39.fasta.gz.unik -o union


    # intersection
    $ unikmer inter Ecoli-MG1655.fasta.gz.unik Ecoli-IAI39.fasta.gz.unik -o inter


    # difference
    $ unikmer diff -t 4 Ecoli-MG1655.fasta.gz.unik Ecoli-IAI39.fasta.gz.unik -o diff


    # -------------------------------------------------------------------------

    $ ls -lh
    -rw-rw-r--. 1 shenwei shenwei 1.6M Aug  9 23:19 Ecoli-IAI39.fasta.gz
    -rw-rw-r--. 1 shenwei shenwei  25M Aug  9 23:29 Ecoli-IAI39.fasta.gz.unik
    -rw-rw-r--. 1 shenwei shenwei 1.4M Aug  9 23:19 Ecoli-MG1655.fasta.gz
    -rw-rw-r--. 1 shenwei shenwei  23M Aug  9 23:29 Ecoli-MG1655.fasta.gz.unik
    -rw-rw-r--. 1 shenwei shenwei  38M Aug  9 23:32 union.unik
    -rw-rw-r--. 1 shenwei shenwei  35M Aug  9 23:33 inter.unik
    -rw-rw-r--. 1 shenwei shenwei  35M Aug  9 23:34 diff.unik

    $ unikmer view Ecoli-MG1655.fasta.gz.unik | wc -l
    9108538
    $ unikmer view Ecoli-IAI39.fasta.gz.unik | wc -l
    9821960
    $ unikmer view union.unik | wc -l
    14402956
    $ unikmer view inter.unik | wc -l
    4527542
    $ unikmer view diff.unik | wc -l
    4580996


## Contributing

We welcome pull requests, bug fixes and issue reports.

## License

[MIT License](https://github.com/shenwei356/unikmer/blob/master/LICENSE)
