# unikmer

unikmer (unique Kmer) is a golang package and a command-line toolkit for
manipulating [Kmers](https://en.wikipedia.org/wiki/K-mer) (k <= 32)
while NOT recording Kmer frequencies.

Every Kmer (k <= 32) is encoded into `uint64`,
and Kmers are stored in builtin `map` in RAM,
no probabilistic data structures are used (I've tested and abandoned them).

## The package

[![GoDoc](https://godoc.org/github.com/shenwei356/unikmer?status.svg)](https://godoc.org/github.com/shenwei356/unikmer)
[![Go Report Card](https://goreportcard.com/badge/github.com/shenwei356/unikmer)](https://goreportcard.com/report/github.com/shenwei356/unikmer)

The unikmer package provides basic manipulations of unique Kmers (NOT including
Kmer frequencies) and provides serialization methods.

### Installation

    go get -u github.com/shenwei356/unikmer

### Benchmark

    $ go test . -bench=Bench* -benchmem
    goos: linux
    goarch: amd64
    pkg: github.com/shenwei356/unikmer
    BenchmarkEncodeK32-4    20000000   90.9 ns/op    0 B/op    0 allocs/op
    BenchmarkDecodeK32-4    20000000    103 ns/op   32 B/op    1 allocs/op
    BenchmarkRevK32-4       20000000   61.9 ns/op    0 B/op    0 allocs/op
    BenchmarkCompK32-4      30000000   53.9 ns/op    0 B/op    0 allocs/op
    BenchmarkRevCompK32-4   10000000    115 ns/op    0 B/op    0 allocs/op

## The toolkit

`unikmer` is a command-line toolkit providing some functions including counting,
format convertion, set operations and searching on unique Kmers.

### Installation

1. Download [binary files](https://github.com/shenwei356/unikmer/releases).

1. Bioconda (not available now)

        conda install unikmer

### Binary file (.unik)

Kmers are saved in gzipped binary file with file extension `.unik`.

Here's the compression rate comparison with `gzip`.
Plain text of Kmers are exported by `unikmer view` and
compressed by `gzip` (default parameters).

    $ ./compression-rates.sh Ecoli-MG1655.fasta.gz > compression-rates.tsv
    $ csvtk -t pretty compression-rates.tsv
    K    plain       .gz        .unik      .gz(%)    .unik(%)
    1    8           64         69         800.000   862.500
    3    256         159        160        62.109    62.500
    5    6144        2070       2267       33.691    36.898
    7    131072      40153      42178      30.634    32.179
    9    2596420     744013     811277     28.655    31.246
    11   35092056    8862922    11083277   25.256    31.583
    13   107877000   22057558   32476722   20.447    30.105
    15   142791328   25348131   43253150   17.752    30.291
    17   163085652   27351927   47477859   16.772    29.112
    19   181600840   27560992   47963931   15.177    26.412
    21   199931204   27939876   34248156   13.975    17.130
    23   218238336   27757535   30860501   12.719    14.141
    25   236543320   28206395   24413804   11.924    10.321
    27   254849280   28019791   27276395   10.995    10.703
    29   273158820   28261557   21435152   10.346    7.847
    31   291473216   28052530   23774695   9.624     8.157

    # plot
    $ csvtk cut -t -f 1,5,6  compression-rates.tsv \
        | csvtk -t gather  -k group  -v value -f 2,3 \
        | csvtk plot line -t -x 1 -y 3 -g 2 \
            --x-min 3 --x-max 40  --y-max 100 --ylab "compression rate (%)" \
            --title Ecoli-MG1655.fasta.gz \
        > compression-rates.tsv.png

![compression-rates.tsv.png](testdata/compression-rates.tsv.png)



### Commands

1. Counting

        count           count Kmer from FASTA/Q sequences
        subset          extract smaller Kmers from binary file

1. Format convertion

        view            read and output binary format to plain text
        dump            convert plain Kmer text to binary format

1. Set operations

        inter           intersection of multiple binary files
        union           union of multiple binary files
        concat          concatenate multiple binary files without removing duplicates
        diff            set difference of multiple binary files

1. Searching

        grep            search Kmer from binary file

1. Misc

        genautocomplete generate shell autocompletion script
        help            Help about any command
        version         print version information and check for update

### Quick Start


    # memusg is for compute time and RAM usage: https://github.com/shenwei356/memusg

    # counting
    $ memusg -t unikmer count -k 31 Ecoli-MG1655.fasta.gz -o Ecoli-MG1655.fasta.gz

    elapsed time: 6.228s
    peak rss: 430.18 MB

    $ ls -lh Ecoli-MG1655.fasta.gz*
    -rw-rw-r--. 1 shenwei shenwei 1.4M Aug  9 23:19 Ecoli-MG1655.fasta.gz
    -rw-rw-r--. 1 shenwei shenwei  23M Aug  9 23:29 Ecoli-MG1655.fasta.gz.k31.unik


    # view
    $ unikmer view Ecoli-MG1655.fasta.gz.k31.unik | head -n 3
    AGCTTTTCATTCTGACTGCAACGGGCAATAT
    GCTTTTCATTCTGACTGCAACGGGCAATATG
    CTTTTCATTCTGACTGCAACGGGCAATATGT

    $ memusg -t unikmer view Ecoli-MG1655.fasta.gz.k31.unik | wc -l

    elapsed time: 2.908s
    peak rss: 19.34 MB

    9108538


    # union
    $ memusg -t unikmer union Ecoli-MG1655.fasta.gz.k31.unik Ecoli-IAI39.fasta.gz.k31.unik -o union

    elapsed time: 10.103s
    peak rss: 773.04 MB

    # intersection
    $ memusg -t unikmer inter Ecoli-MG1655.fasta.gz.k31.unik Ecoli-IAI39.fasta.gz.k31.unik -o inter

    elapsed time: 7.955s
    peak rss: 400.71 MB


    # difference
    $ memusg -t unikmer diff -t 1 Ecoli-MG1655.fasta.gz.k31.unik Ecoli-IAI39.fasta.gz.k31.unik -o diff

    elapsed time: 8.137s
    peak rss: 400.45 MB

    # -------------------------------------------------------------------------

    $ ls -lh
    -rw-rw-r--. 1 shenwei shenwei 1.6M Aug  9 23:19 Ecoli-IAI39.fasta.gz
    -rw-rw-r--. 1 shenwei shenwei  25M Aug  9 23:29 Ecoli-IAI39.fasta.gz.k31.unik
    -rw-rw-r--. 1 shenwei shenwei 1.4M Aug  9 23:19 Ecoli-MG1655.fasta.gz
    -rw-rw-r--. 1 shenwei shenwei  23M Aug  9 23:29 Ecoli-MG1655.fasta.gz.k31.unik
    -rw-rw-r--. 1 shenwei shenwei  38M Aug  9 23:32 union.k31.unik
    -rw-rw-r--. 1 shenwei shenwei  35M Aug  9 23:33 inter.k31.unik
    -rw-rw-r--. 1 shenwei shenwei  35M Aug  9 23:34 diff.k31.unik

    $ unikmer view Ecoli-MG1655.fasta.gz.k31.unik | wc -l
    9108538
    $ unikmer view Ecoli-IAI39.fasta.gz.k31.unik | wc -l
    9821960
    $ unikmer view union.k31.unik | wc -l
    14402956
    $ unikmer view inter.k31.unik | wc -l
    4527542
    $ unikmer view diff.k31.unik | wc -l
    4580996


## Contributing

We welcome pull requests, bug fixes and issue reports.

## License

[MIT License](https://github.com/shenwei356/unikmer/blob/master/LICENSE)
