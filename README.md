# unikmer: Toolkit for k-mer with taxonomic information

`unikmer` is a toolkit for nucleic acid [k-mer](https://en.wikipedia.org/wiki/K-mer) analysis, 
providing functions
including set operation k-mers (sketch) optional with
TaxIds but without count information.

K-mers are either encoded (k<=32) or hashed (arbitrary k) into `uint64`,
and serialized in binary file with extension `.unik`.

TaxIds can be assigned when counting k-mers from genome sequences,
and LCA (Lowest Common Ancestor) is computed during set opertions
including computing union, intersecton, set difference, unique and
repeated k-mers.

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
## Table of Contents

- [Using cases](#using-cases)
- [Installation](#installation)
- [Commands](#commands)
- [Binary file](#binary-file)
- [Quick start](#quick-start)
- [Support](#support)
- [License](#license)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

## Using cases

- Finding conserved regions in all genomes of a species.
- Finding species/strain-specific sequences for designing probes/primers.

## Installation

1. Downloading [executable binary files](https://github.com/shenwei356/unikmer/releases).

1. Via Bioconda [![Anaconda Cloud](https://anaconda.org/bioconda/unikmer/badges/version.svg)](https://anaconda.org/bioconda/unikmer) [![downloads](https://anaconda.org/bioconda/unikmer/badges/downloads.svg)](https://anaconda.org/bioconda/unikmer)

        conda install -c bioconda unikmer

## Commands

[Usages](https://bioinf.shenwei.me/unikmer/usage)

1. Counting

        count           Generate k-mers (sketch) from FASTA/Q sequences

1. Information

        info            Information of binary files
        num             Quickly inspect number of k-mers in binary files

1. Format conversion

        view            Read and output binary format to plain text
        dump            Convert plain k-mer text to binary format

        encode          Encode plain k-mer text to integer
        decode          Decode encoded integer to k-mer text
        

1. Set operations

        concat          Concatenate multiple binary files without removing duplicates
        inter           Intersection of multiple binary files
        common          Find k-mers shared by most of multiple binary files
        union           Union of multiple binary files
        diff            Set difference of multiple binary files

1. Split and merge

        sort            Sort k-mers in binary files to reduce file size
        split           Split k-mers into sorted chunk files
        tsplit          Split k-mers according to TaxId
        merge           Merge k-mers from sorted chunk files

1. Subset

        head            Extract the first N k-mers
        sample          Sample k-mers from binary files
        grep            Search k-mers from binary files
        filter          Filter out low-complexity k-mers
        rfilter         Filter k-mers by taxonomic rank

1. Searching on genomes

        locate          Locate k-mers in genome
        uniqs           Mapping k-mers back to genome and find unique subsequences

1. Misc

        autocompletion  Generate shell autocompletion script
        version         Print version information and check for update

## Binary file

[![Go Reference](https://pkg.go.dev/badge/github.com/shenwei356/unik.svg)](https://pkg.go.dev/github.com/shenwei356/unik)

K-mers (represented in `uint64` in RAM ) are serialized in 8-Byte
(or less Bytes for shorter k-mers in compact format,
or much less Bytes for sorted k-mers) arrays and
optionally compressed in gzip format with extension of `.unik`.
TaxIds are optionally stored next to k-mers with 4 or less bytes.

### Compression ratio comparison

No TaxIds stored in this test.

![cr.jpg](testdata/cr.jpg)

label           |encoded-kmer<sup>a</sup>|gzip-compressed<sup>b</sup>|compact-format<sup>c</sup>|sorted<sup>d</sup>|comment
:---------------|:----------------------:|:-------------------------:|:------------------------:|:----------------:|:------------------------------------------------------
`plain`         |                        |                           |                          |                  |plain text
`gzip`          |                        |✔                          |                          |                  |gzipped plain text
`unik.default`  |✔                       |✔                          |                          |                  |gzipped encoded k-mers in fixed-length byte array
`unik.compat`   |✔                       |✔                          |✔                         |                  |gzipped encoded k-mers in shorter fixed-length byte array
`unik.sorted`   |✔                       |✔                          |                          |✔                 |gzipped sorted encoded k-mers


- <sup>a</sup> One k-mer is encoded as `uint64` and serialized in 8 Bytes.
- <sup>b</sup> K-mers file is compressed in gzip format by default,
  users can switch on global option `-C/--no-compress` to output non-compressed file.
- <sup>c</sup> One k-mer is encoded as `uint64` and serialized in 8 Bytes by default.
 However few Bytes are needed for short k-mers, e.g., 4 Bytes are enough for
  15-mers (30 bits). This makes the file more compact with smaller file size,
  controled by global option `-c/--compact `.
- <sup>d</sup> One k-mer is encoded as `uint64`, all k-mers are sorted and compressed
  using varint-GB algorithm.
- In all test, flag `--canonical` is ON when running `unikmer count`.


## Quick Start


    # memusg is for compute time and RAM usage: https://github.com/shenwei356/memusg


    # counting (only keep the canonical k-mers and compact output)
    # memusg -t unikmer count -k 23 Ecoli-IAI39.fasta.gz -o Ecoli-IAI39.fasta.gz.k23 --canonical --compact
    $ memusg -t unikmer count -k 23 Ecoli-MG1655.fasta.gz -o Ecoli-MG1655.fasta.gz.k23 --canonical --compact
    elapsed time: 0.897s
    peak rss: 192.41 MB


    # counting (only keep the canonical k-mers and sort k-mers)
    # memusg -t unikmer count -k 23 Ecoli-IAI39.fasta.gz -o Ecoli-IAI39.fasta.gz.k23.sorted --canonical --sort
    $ memusg -t unikmer count -k 23 Ecoli-MG1655.fasta.gz -o Ecoli-MG1655.fasta.gz.k23.sorted --canonical --sort
    elapsed time: 1.136s
    peak rss: 227.28 MB
    
    
    # counting and assigning global TaxIds
    $ unikmer count -k 23 -K -s Ecoli-IAI39.fasta.gz -o Ecoli-IAI39.fasta.gz.k23.sorted   -t 585057
    $ unikmer count -k 23 -K -s Ecoli-MG1655.fasta.gz -o Ecoli-MG1655.fasta.gz.k23.sorted -t 511145
    $ unikmer count -k 23 -K -s A.muciniphila-ATCC_BAA-835.fasta.gz -o A.muciniphila-ATCC_BAA-835.fasta.gz.sorted -t 349741
    
    # counting minimizer and ouputting in linear order
    $ unikmer count -k 23 -W 5 -H -K -l A.muciniphila-ATCC_BAA-835.fasta.gz -o A.muciniphila-ATCC_BAA-835.fasta.gz.m

    # view
    $ unikmer view Ecoli-MG1655.fasta.gz.k23.sorted.unik --show-taxid | head -n 3
    AAAAAAAAACCATCCAAATCTGG 511145
    AAAAAAAAACCGCTAGTATATTC 511145
    AAAAAAAAACCTGAAAAAAACGG 511145
    
    # view (hashed k-mers needs original FASTA/Q file)
    $ unikmer view --show-code --genome A.muciniphila-ATCC_BAA-835.fasta.gz A.muciniphila-ATCC_BAA-835.fasta.gz.m.unik | head -n 3
    CATCCGCCATCTTTGGGGTGTCG 1210726578792
    AGCGCAAAATCCCCAAACATGTA 2286899379883
    AACTGATTTTTGATGATGACTCC 3542156397282
    
    # find the positions of k-mers
    $ seqkit locate -M A.muciniphila-ATCC_BAA-835.fasta.gz \
        -f <(unikmer view -a -g A.muciniphila-ATCC_BAA-835.fasta.gz A.muciniphila-ATCC_BAA-835.fasta.gz.m.unik | seqkit head -n 5 ) \
        | csvtk sort -t -k start:n | head -n 6 | csvtk pretty -t
    seqID         patternName           pattern                   strand   start   end
    -----------   -------------------   -----------------------   ------   -----   ---
    NC_010655.1   2090893901864583115   ATCTTATAAAATAACCACATAAC   +        3       25
    NC_010655.1   696051979077366638    TTATAAAATAACCACATAACTTA   +        6       28
    NC_010655.1   390297872016815006    TATAAAATAACCACATAACTTAA   +        7       29
    NC_010655.1   2582400417208090837   AAAATAACCACATAACTTAAAAA   +        10      32
    NC_010655.1   3048591415312050785   TAACCACATAACTTAAAAAGAAT   +        14      36
    
    # stats
    $ unikmer stats *.unik -a -j 10
    file                                              k  canonical  hashed  scaled  include-taxid  global-taxid  sorted  compact  gzipped  version     number  description
    A.muciniphila-ATCC_BAA-835.fasta.gz.m.unik       23          ✓       ✓       ✕              ✕                     ✕        ✕        ✓     v5.0    860,900  
    A.muciniphila-ATCC_BAA-835.fasta.gz.sorted.unik  23          ✓       ✕       ✕              ✕        349741       ✓        ✕        ✓     v5.0  2,630,905  
    Ecoli-IAI39.fasta.gz.k23.sorted.unik             23          ✓       ✕       ✕              ✕        585057       ✓        ✕        ✓     v5.0  4,902,266  
    Ecoli-IAI39.fasta.gz.k23.unik                    23          ✓       ✕       ✕              ✕                     ✕        ✓        ✓     v5.0  4,902,266  
    Ecoli-MG1655.fasta.gz.k23.sorted.unik            23          ✓       ✕       ✕              ✕        511145       ✓        ✕        ✓     v5.0  4,546,632  
    Ecoli-MG1655.fasta.gz.k23.unik                   23          ✓       ✕       ✕              ✕                     ✕        ✓        ✓     v5.0  4,546,632 

    
    # concat
    $ memusg -t unikmer concat *.k23.sorted.unik -o concat.k23 -c
    elapsed time: 1.020s
    peak rss: 25.86 MB


    
    # union
    $ memusg -t unikmer union *.k23.sorted.unik -o union.k23 -s
    elapsed time: 3.991s
    peak rss: 590.92 MB
    
    
    # or sorting with limited memory.
    # note that taxonomy database need some memory.
    $ memusg -t unikmer sort *.k23.sorted.unik -o union2.k23 -u -m 1M
    elapsed time: 3.538s
    peak rss: 324.2 MB
    
    $ unikmer view -t union.k23.unik | md5sum 
    4c038832209278840d4d75944b29219c  -
    $ unikmer view -t union2.k23.unik | md5sum 
    4c038832209278840d4d75944b29219c  -
    
    
    # duplicate k-mers
    $ memusg -t unikmer sort *.k23.sorted.unik -o dup.k23 -d -m 1M
    elapsed time: 1.143s
    peak rss: 240.18 MB

    
    # intersection
    $ memusg -t unikmer inter *.k23.sorted.unik -o inter.k23
    elapsed time: 1.481s
    peak rss: 399.94 MB
    

    # difference
    $ memusg -t unikmer diff -j 10 *.k23.sorted.unik -o diff.k23 -s
    elapsed time: 0.793s
    peak rss: 338.06 MB


    $ ls -lh *.unik
    -rw-r--r-- 1 shenwei shenwei 9.5M  2月 13 00:55 A.muciniphila-ATCC_BAA-835.fasta.gz.sorted.unik
    -rw-r--r-- 1 shenwei shenwei  46M  2月 13 00:59 concat.k23.unik
    -rw-r--r-- 1 shenwei shenwei 8.7M  2月 13 01:00 diff.k23.unik
    -rw-r--r-- 1 shenwei shenwei  11M  2月 13 01:04 dup.k23.unik
    -rw-r--r-- 1 shenwei shenwei  18M  2月 13 00:55 Ecoli-IAI39.fasta.gz.k23.sorted.unik
    -rw-r--r-- 1 shenwei shenwei  21M  2月 13 00:48 Ecoli-IAI39.fasta.gz.k23.unik
    -rw-r--r-- 1 shenwei shenwei  17M  2月 13 00:55 Ecoli-MG1655.fasta.gz.k23.sorted.unik
    -rw-r--r-- 1 shenwei shenwei  19M  2月 13 00:48 Ecoli-MG1655.fasta.gz.k23.unik
    -rw-r--r-- 1 shenwei shenwei 9.5M  2月 13 00:59 inter.k23.unik
    -rw-r--r-- 1 shenwei shenwei  27M  2月 13 01:04 union2.k23.unik
    -rw-r--r-- 1 shenwei shenwei  27M  2月 13 00:58 union.k23.unik


    $ unikmer stats *.unik -a -j 10
    file                                              k  canonical  hashed  scaled  include-taxid  global-taxid  sorted  compact  gzipped  version     number  description
    A.muciniphila-ATCC_BAA-835.fasta.gz.m.unik       23          ✓       ✓       ✕              ✕                     ✕        ✕        ✓     v5.0    860,900  
    A.muciniphila-ATCC_BAA-835.fasta.gz.sorted.unik  23          ✓       ✕       ✕              ✕        349741       ✓        ✕        ✓     v5.0  2,630,905  
    concat.k23.unik                                  23          ✓       ✕       ✕              ✓                     ✕        ✓        ✓     v5.0         -1  
    diff.k23.unik                                    23          ✓       ✕       ✕              ✓                     ✕        ✕        ✓     v5.0  2,326,096  
    dup.k23.unik                                     23          ✓       ✕       ✕              ✓                     ✓        ✕        ✓     v5.0          0  
    Ecoli-IAI39.fasta.gz.k23.sorted.unik             23          ✓       ✕       ✕              ✕        585057       ✓        ✕        ✓     v5.0  4,902,266  
    Ecoli-IAI39.fasta.gz.k23.unik                    23          ✓       ✕       ✕              ✕                     ✕        ✓        ✓     v5.0  4,902,266  
    Ecoli-MG1655.fasta.gz.k23.sorted.unik            23          ✓       ✕       ✕              ✕        511145       ✓        ✕        ✓     v5.0  4,546,632  
    Ecoli-MG1655.fasta.gz.k23.unik                   23          ✓       ✕       ✕              ✕                     ✕        ✓        ✓     v5.0  4,546,632  
    inter.k23.unik                                   23          ✓       ✕       ✕              ✓                     ✓        ✕        ✓     v5.0  2,576,170  
    union2.k23.unik                                  23          ✓       ✕       ✕              ✓                     ✓        ✕        ✓     v5.0  6,872,728  
    union.k23.unik                                   23          ✓       ✕       ✕              ✓                     ✓        ✕        ✓     v5.0  6,872,728


    # -----------------------------------------------------------------------------------------

    # mapping k-mers to genome
    g=Ecoli-IAI39.fasta
    f=inter.k23.unik

    # to fasta
    unikmer view $f -a -o $f.fa.gz

    # make index
    bwa index $g; samtools faidx $g

    ncpu=12
    ls $f.fa.gz \
        | rush -j 1 -v ref=$g -v j=$ncpu \
            'bwa aln -o 0 -l 17 -k 0 -t {j} {ref} {} \
                | bwa samse {ref} - {} \
                | samtools view -bS > {}.bam; \
             samtools sort -T {}.tmp -@ {j} {}.bam -o {}.sorted.bam; \
             samtools index {}.sorted.bam; \
             samtools flagstat {}.sorted.bam > {}.sorted.bam.flagstat; \
             /bin/rm {}.bam '  

## Support

Please [open an issue](https://github.com/shenwei356/unikmer/issues) to report bugs,
propose new functions or ask for help.

## License

[MIT License](https://github.com/shenwei356/unikmer/blob/master/LICENSE)
