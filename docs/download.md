# Download

unikmer is implemented in [Go](https://golang.org/) programming language,
statically-linked executable binary files are [freely available](https://github.com/shenwei356/unikmer/releases).

## Current Version

### [v0.19.1](https://github.com/shenwei356/unikmer/releases/tag/v0.19.1) - 2022-12-26 [![Github Releases (by Release)](https://img.shields.io/github/downloads/shenwei356/unikmer/v0.19.1/total.svg)](https://github.com/shenwei356/unikmer/releases/tag/v0.19.1)

- `unikmer`: When environment variable `UNIKMER_DB` is set, explicitly setting `--data-dir` will overide the value of `UNIKMER_DB`.
- `unikmer uniqs`: skip sequences shorter than K.
- `unikmer count/encode`: [limit the maximum k-mer size to 64](https://github.com/bcgsc/ntHash/issues/41).

### Links

OS     |Arch      |File, 中国镜像                                                                                                                                                                                  |Download Count
:------|:---------|:------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|:-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Linux  |**64-bit**|[**unikmer_linux_amd64.tar.gz**](https://github.com/shenwei356/unikmer/releases/download/v0.19.1/unikmer_linux_amd64.tar.gz), <br/> [中国镜像](http://app.shenwei.me/data/unikmer/unikmer_linux_amd64.tar.gz)                  |[![Github Releases (by Asset)](https://img.shields.io/github/downloads/shenwei356/unikmer/latest/unikmer_linux_amd64.tar.gz.svg?maxAge=3600)](https://github.com/shenwei356/unikmer/releases/download/v0.19.1/unikmer_linux_amd64.tar.gz)
Linux  |arm64|[**unikmer_linux_arm64.tar.gz**](https://github.com/shenwei356/unikmer/releases/download/v0.19.1/unikmer_linux_arm64.tar.gz), <br/> [中国镜像](http://app.shenwei.me/data/unikmer/unikmer_linux_arm64.tar.gz)                  |[![Github Releases (by Asset)](https://img.shields.io/github/downloads/shenwei356/unikmer/latest/unikmer_linux_arm64.tar.gz.svg?maxAge=3600)](https://github.com/shenwei356/unikmer/releases/download/v0.19.1/unikmer_linux_arm64.tar.gz)
macOS  |**64-bit**|[**unikmer_darwin_amd64.tar.gz**](https://github.com/shenwei356/unikmer/releases/download/v0.19.1/unikmer_darwin_amd64.tar.gz), <br/> [中国镜像](http://app.shenwei.me/data/unikmer/unikmer_darwin_amd64.tar.gz)               |[![Github Releases (by Asset)](https://img.shields.io/github/downloads/shenwei356/unikmer/latest/unikmer_darwin_amd64.tar.gz.svg?maxAge=3600)](https://github.com/shenwei356/unikmer/releases/download/v0.19.1/unikmer_darwin_amd64.tar.gz)
macOS  |arm64|[**unikmer_darwin_arm64.tar.gz**](https://github.com/shenwei356/unikmer/releases/download/v0.19.1/unikmer_darwin_arm64.tar.gz), <br/> [中国镜像](http://app.shenwei.me/data/unikmer/unikmer_darwin_arm64.tar.gz)               |[![Github Releases (by Asset)](https://img.shields.io/github/downloads/shenwei356/unikmer/latest/unikmer_darwin_arm64.tar.gz.svg?maxAge=3600)](https://github.com/shenwei356/unikmer/releases/download/v0.19.1/unikmer_darwin_arm64.tar.gz)
Windows|**64-bit**|[**unikmer_windows_amd64.exe.tar.gz**](https://github.com/shenwei356/unikmer/releases/download/v0.19.1/unikmer_windows_amd64.exe.tar.gz), <br/> [中国镜像](http://app.shenwei.me/data/unikmer/unikmer_windows_amd64.exe.tar.gz)|[![Github Releases (by Asset)](https://img.shields.io/github/downloads/shenwei356/unikmer/latest/unikmer_windows_amd64.exe.tar.gz.svg?maxAge=3600)](https://github.com/shenwei356/unikmer/releases/download/v0.19.1/unikmer_windows_amd64.exe.tar.gz)

*Notes:*

- please open an issue to request binaries for other platforms.
- run `unikmer version` to check update !!!
- run `unikmer autocompletion` to update shell autocompletion script !!!


## Installation 

#### Method 1: Install using conda [![Anaconda Cloud](https://anaconda.org/bioconda/unikmer/badges/version.svg)](https://anaconda.org/bioconda/unikmer) [![downloads](https://anaconda.org/bioconda/unikmer/badges/downloads.svg)](https://anaconda.org/bioconda/unikmer)

    conda install -c bioconda unikmer

#### Method 2: Download binaries

[Download](https://github.com/shenwei356/unikmer/releases) the compressed
executable file of your operating system,
and decompress it with `tar -zxvf *.tar.gz` command or other tools.
And then:

- **For Linux-like systems**
    - If you have root privilege, simply copy it to `/usr/local/bin`:

            sudo cp unikmer /usr/local/bin/

    - Or copy to anywhere in the environment variable `PATH`:

            mkdir -p $HOME/bin/; cp unikmer $HOME/bin/

- **For Windows**, just copy `unikmer.exe` to `C:\WINDOWS\system32`.

#### Method 3: Compile from source

1. [Install go](https://go.dev/doc/install)

        wget https://go.dev/dl/go1.17.13.linux-amd64.tar.gz

        tar -zxf go1.17.13.linux-amd64.tar.gz -C $HOME/

        # or
        #   echo "export PATH=$PATH:$HOME/go/bin" >> ~/.bashrc
        #   source ~/.bashrc
        export PATH=$PATH:$HOME/go/bin

2. Compile KMCP

        # ------------- the latest stable version -------------

        go get -v -u github.com/shenwei356/unikmer/unikmer

        # The executable binary file is located in:
        #   ~/go/bin/unikmer
        # You can also move it to anywhere in the $PATH
        mkdir -p $HOME/bin
        cp ~/go/bin/unikmer $HOME/bin/


        # --------------- the development version --------------

        git clone https://github.com/shenwei356/unikmer
        cd unikmer/unikmer/
        go build

        # The executable binary file is located in:
        #   ./unikmer
        # You can also move it to anywhere in the $PATH
        mkdir -p $HOME/bin
        cp ./unikmer $HOME/bin/


## Shell-completion

Supported shell: bash|zsh|fish|powershell

Bash:

    # generate completion shell
    unikmer autocompletion --shell bash

    # configure if never did.
    # install bash-completion if the "complete" command is not found.
    echo "for bcfile in ~/.bash_completion.d/* ; do source \$bcfile; done" >> ~/.bash_completion
    echo "source ~/.bash_completion" >> ~/.bashrc

Zsh:

    # generate completion shell
    unikmer autocompletion --shell zsh --file ~/.zfunc/_unikmer

    # configure if never did
    echo 'fpath=( ~/.zfunc "${fpath[@]}" )' >> ~/.zshrc
    echo "autoload -U compinit; compinit" >> ~/.zshrc

fish:

    unikmer autocompletion --shell fish --file ~/.config/fish/completions/unikmer.fish

## Release History

### [v0.19.0](https://github.com/shenwei356/unikmer/releases/tag/v0.19.0) - 2022-04-25 [![Github Releases (by Release)](https://img.shields.io/github/downloads/shenwei356/unikmer/v0.19.0/total.svg)](https://github.com/shenwei356/unikmer/releases/tag/v0.19.0)

- rename command `genautocomplete` to `autocompletion`.
- remove command `help`.
- change default value of option `-j` from `2` to `4`.
- `unikmer count/uniqs/locate`: new flag `-B/--seq-name-filter` for filtering out unwanted sequences like plasmids.
- `unikmer count`: add support of `.xz` and `.zst` files.

### [v0.18.8](https://github.com/shenwei356/unikmer/releases/tag/v0.18.8) - 2021-09-17

- use new version of nthash with better performance.
- `unikmer info`: fix typoes.

### [v0.18.7](https://github.com/shenwei356/unikmer/releases/tag/v0.18.7) - 2021-08-30 [![Github Releases (by Release)](https://img.shields.io/github/downloads/shenwei356/unikmer/v0.18.7/total.svg)](https://github.com/shenwei356/unikmer/releases/tag/v0.18.7)

- `unikmer`: better counting speed by upstream optimization of FASTA/Q parsing.
- `unikmer concat`: fix parsing flag `-n`.
  
### [v0.17.3](https://github.com/shenwei356/unikmer/releases/tag/v0.17.3) - 2021-05-16

- `unikmer`: fix buiding for 386. #21

### [v0.17.2](https://github.com/shenwei356/unikmer/releases/tag/v0.17.2) - 2021-02-05 [![Github Releases (by Release)](https://img.shields.io/github/downloads/shenwei356/unikmer/v0.17.2/total.svg)](https://github.com/shenwei356/unikmer/releases/tag/v0.17.2)

- `unikmer`: slightly speedup for computing LCA.
- `unikmer rfilter:` 
    - flag `-E/--equal-to` supports multiple values.
    - new flag `-n/--save-predictable-norank`: do not discard some special ranks without order when using -L, where rank of the closest higher node is still lower than rank cutoff.

### [v0.17.1](https://github.com/shenwei356/unikmer/releases/tag/v0.17.1) - 2021-01-18

- `unikmer rfilter:` change handling of black list.

### [v0.17.0](https://github.com/shenwei356/unikmer/releases/tag/v0.17.0) - 2021-01-15

- **syncmer value changed with different hash method**.
- `unikmer count`: syncmer value changed.

### [v0.16.1](https://github.com/shenwei356/unikmer/releases/tag/v0.16.1) - 2020-12-28 [![Github Releases (by Release)](https://img.shields.io/github/downloads/shenwei356/unikmer/v0.16.1/total.svg)](https://github.com/shenwei356/unikmer/releases/tag/v0.16.1)

- change Header.Number from `int64` to `uint64`
- `unikmer info`: fix recounting problem for unsorted kmers but with Number.
  
### [v0.16.0](https://github.com/shenwei356/unikmer/releases/tag/v0.16.0) - 2020-12-28

- `unikmer`:
    - **binary file format change**: fix reading long description, and bump version to `5.0`.
    - better binary file parsing performance.

### [v0.15.0](https://github.com/shenwei356/unikmer/releases/tag/v0.15.0) - 2020-12-25 [![Github Releases (by Release)](https://img.shields.io/github/downloads/shenwei356/unikmer/v0.15.0/total.svg)](https://github.com/shenwei356/unikmer/releases/tag/v0.15.0)

- `unikmer`:
    - binary file minor change: increase description maximal length from 128 B to 1KB.
    - separating k-mers (sketches) indexing and searching from `unikmer`, including `unikmer db info/index/search`.
- `unikmer count`: fix syncmer.
- `unikmer dump`: new flag `--hashed`.
- rename `unikmer stats` to `unikmer info`, and add new column `description`.
  
### [v0.14.0](https://github.com/shenwei356/unikmer/releases/tag/v0.14.0) - 2020-11-25 [![Github Releases (by Release)](https://img.shields.io/github/downloads/shenwei356/unikmer/v0.14.0/total.svg)](https://github.com/shenwei356/unikmer/releases/tag/v0.14.0)

- `unikmer union`: fix bug when flag `-s` not given.
- `unikmer count/uniqs/locate`: performance improvement on generating k-mers.
- `unikmer count/db`: support scaled/minizimer/syncmer sketch.
- `unikmer stats`: change format.
  
### [v0.13.0](https://github.com/shenwei356/unikmer/releases/tag/v0.13.0) - 2020-10-23 [![Github Releases (by Release)](https://img.shields.io/github/downloads/shenwei356/unikmer/v0.13.0/total.svg)](https://github.com/shenwei356/unikmer/releases/tag/v0.13.0)

- new command `unikmer common`: Finding k-mers shared by most of multiple binary files.
- `unikmer common/count/diff/grep/rfilter/sort/split/union`: faster sorting.
- `unikmer uniqs`: better result for flag `--circular`.
- `unikmer search`: fix a bug when searching on database with more than one hash.
  
### [v0.12.0](https://github.com/shenwei356/unikmer/releases/tag/v0.12.0) - 2020-09-24 [![Github Releases (by Release)](https://img.shields.io/github/downloads/shenwei356/unikmer/v0.12.0/total.svg)](https://github.com/shenwei356/unikmer/releases/tag/v0.12.0)

- `unikmer`:
    - support longer k (k>32) by saving ntHash.
    - new flag `-nocheck-file` for not checking binary file.
- new commands:
    - `unikmer db index`: constructing index from binary files
    - `unikmer db info`: printing information of index file
    - `unikmer db search`: searching sequence from index database
- `unikmer rfilter`: change format of rank order file.
- `unikmer inter/union`: speedup for single input file.
- `unikmer concat`:
    - new flag `-t/--taxid` for assigning global taxid, this can slightly reduce file size.
    - new flag `-n/--number` for setting number of k-mers.
- `unikmer num`:
    - new flag `-f/--force` for counting k-mers.
- `unikmer locate`: output in BED6.
- `unikmer locate/uniqs`: support multiple genome files.
- `unikmer uniqs`:
    - stricter multiple mapping limit.
    - new flag `-W/--seqs-in-a-file-as-one-genome`.
- `unikmer count`:
    - new flag `-u/--unique` for output unique (single copy) kmers

### [v0.11.0](https://github.com/shenwei356/unikmer/releases/tag/v0.11.0) - 2020-07-06 [![Github Releases (by Release)](https://img.shields.io/github/downloads/shenwei356/unikmer/v0.11.0/total.svg)](https://github.com/shenwei356/unikmer/releases/tag/v0.11.0)

- new command: `unikmer rfilter` for filtering k-mers by taxonomic rank.
- `unikmer inter`: new flag `-m/--mix-taxid` allowing part of files being whithout taxids.
- `unikmer dump`: fix a nil pointer bug.
- `unikmer count`:
    - fix checking taxid in sequence header.
    - fix setting global taxid.
- `unikmer count/diff/union`: slightly reduce memory and speedup when sorting k-mers.
- `unikmer filter`: change scoring.
- `unikmer count/locate/uniqs`: remove flag `--circular`.

### [v0.10.0](https://github.com/shenwei356/unikmer/releases/tag/v0.10.0) - 2020-05-21 [![Github Releases (by Release)](https://img.shields.io/github/downloads/shenwei356/unikmer/v0.10.0/total.svg)](https://github.com/shenwei356/unikmer/releases/tag/v0.10.0)

- `unikmer`: fix loading custom taxonomy files.
- `unikmer count`:
    - new flag `-d` for only count duplicate k-mers, for removing singleton in FASTQ.
    - fix nil pointer bug of `-t`.
- `unikmer split`: fix memery and last odd k-mer mising bug for given ONE sorted input file.
- `unikmer sort`: skip loading taxonomy data when neither `-u` or `-d` given.
- `unikmer diff`: 2X speedup, and requiring 1th file being sorted.
- `unikmer inter`: 2-5X speedup, and requiring all files being sorted, sorted output by default.

### [v0.9.0](https://github.com/shenwei356/unikmer/releases/tag/v0.9.0) - 2020-02-18 [![Github Releases (by Release)](https://img.shields.io/github/downloads/shenwei356/unikmer/v0.9.0/total.svg)](https://github.com/shenwei356/unikmer/releases/tag/v0.9.0)
 
- `unikmer`: **new binary format supporting optional Taxids**.
- deleted command: `unikmer subset`.
- new command: `unikmer head` for extracting the first N k-mers.
- new command: `unikmer tsplit` for splitting k-mers according to taxid.
- `unikmer grep`: support searching with taxids.
- `unikmer count`: support parsing taxid from FASTA/Q header.
  
### [v0.8.0](https://github.com/shenwei356/unikmer/releases/tag/v0.8.0) - 2019-02-09 [![Github Releases (by Release)](https://img.shields.io/github/downloads/shenwei356/unikmer/v0.8.0/total.svg)](https://github.com/shenwei356/unikmer/releases/tag/v0.8.0)

- `unikmer`:
    - new option `-i/--infile-list`, if given, files in the list file are appended to files from cli arguments.
    - improve performance of binary file reading and writing.
- `unikmer sort/split/merge`: safer forcing deletion of existed outdir, and better log.
- `unikmer split`: performance improvement for single sorted input file.
- `unikmer sort`: performance improvement for using `-m/--chunk-size`.
- `unikmer grep`: rewrite, support loading queries from .unik files.
- `unikmer dump`: fix number information in output file.
- `unikmer concat`: new flag `-s/--sorted`.
  
### [v0.7.0](https://github.com/shenwei356/unikmer/releases/tag/v0.7.0) - 2019-09-29 [![Github Releases (by Release)](https://img.shields.io/github/downloads/shenwei356/unikmer/v0.7.0/total.svg)](https://github.com/shenwei356/unikmer/releases/tag/v0.7.0)

- new command `unikmer filter`: filter low-complexity k-mers.
- new command `unikmer split`: split k-mers into sorted chunk files.
- new command `unikmer merge`: merge from sorted chunk files.
- `unikmer view`:
    - new option `-N/--show-code-only` for only showing encoded integers.
    - fix output error for `-q/--fastq`.
- `unikmer uniqs`:
    - new option `-x/--max-cont-non-uniq-kmers` for limiting max continuous non-unique k-mers.
    - new option `-X/--max-num-cont-non-uniq-kmers` for limiting max number of continuous non-unique k-mers.
    - fix bug for `-m/--min-len`.
- `unikmer union`:
    - new option `-d/--repeated` for only printing duplicate k-mers.
- `unikmer sort`:
    - new option `-u/--unique` for removing duplicate k-mers.
    - new option `-d/--repeated` for only printing duplicate k-mers.
    - new option `-m/--chunk-size` for limiting maximum memory for sorting.
- `unikmer diff`:
    - small speed improvements.

### [v0.6.2](https://github.com/shenwei356/unikmer/releases/tag/v0.6.2) - 2019-01-21 [![Github Releases (by Release)](https://img.shields.io/github/downloads/shenwei356/unikmer/v0.6.2/total.svg)](https://github.com/shenwei356/unikmer/releases/tag/v0.6.2)

- `unikmer encode`: better output for bits presentation of encoded k-mers (`-a/--all`)

### [v0.6.1](https://github.com/shenwei356/unikmer/releases/tag/v0.6.1) - 2019-01-21 [![Github Releases (by Release)](https://img.shields.io/github/downloads/shenwei356/unikmer/v0.6.1/total.svg)](https://github.com/shenwei356/unikmer/releases/tag/v0.6.1)

- `unikmer dump`: 
    - new option `-K/--canonical` to keep the canonical k-mers.
    - new option `-k/--canonical-only` to only keep the canonical k-mers.
    - new option `-s/--sorted` to save sorted k-mers.
- `unikmer encode`: add option `-K/--canonical` to keep the canonical k-mers.
  
### [v0.6.0](https://github.com/shenwei356/unikmer/releases/tag/v0.6.0) - 2019-01-20 [![Github Releases (by Release)](https://img.shields.io/github/downloads/shenwei356/unikmer/v0.6.0/total.svg)](https://github.com/shenwei356/unikmer/releases/tag/v0.6.0)

- `unikmer`: check encoded integer overflow
- new command `unikmer encode`: encode plain k-mer text to integer
- new command `unikmer decode`: decode encoded integer to k-mer text
  
### [v0.5.3](https://github.com/shenwei356/unikmer/releases/tag/v0.5.3) - 2018-11-28

- `unikmer count/dump`: check file before handling them.

### [v0.5.2](https://github.com/shenwei356/unikmer/releases/tag/v0.5.2) - 2018-11-28 [![Github Releases (by Release)](https://img.shields.io/github/downloads/shenwei356/unikmer/v0.5.2/total.svg)](https://github.com/shenwei356/unikmer/releases/tag/v0.5.2)

- `unikmer locate`: fix bug.
- `unikmer`: doc update.

### [v0.5.1](https://github.com/shenwei356/unikmer/releases/tag/v0.5.1) - 2018-11-07 [![Github Releases (by Release)](https://img.shields.io/github/downloads/shenwei356/unikmer/v0.5.1/total.svg)](https://github.com/shenwei356/unikmer/releases/tag/v0.5.1)

- `unikmer locate/uniqs`: fix options checking.

### [v0.5.0](https://github.com/shenwei356/unikmer/releases/tag/v0.5.0) - 2018-11-07

- `unikmer diff`: fix concurrency bug when cloning kmers from first file.
- new command `unikmer locate`: locate Kmers in genome.
- new command `unikmer uniqs`: mapping Kmers back to genome and find unique subsequences.
  
### [v0.4.4](https://github.com/shenwei356/unikmer/releases/tag/v0.4.4) - 2018-10-27 [![Github Releases (by Release)](https://img.shields.io/github/downloads/shenwei356/unikmer/v0.4.4/total.svg)](https://github.com/shenwei356/unikmer/releases/tag/v0.4.4)

- `unikmer`: add global option `-L/--compression-level`.
- `unikmer diff`: reduce memory occupation, speed not affected.

### [v0.4.3](https://github.com/shenwei356/unikmer/releases/tag/v0.4.3) - 2018-10-13 [![Github Releases (by Release)](https://img.shields.io/github/downloads/shenwei356/unikmer/v0.4.3/total.svg)](https://github.com/shenwei356/unikmer/releases/tag/v0.4.3)

- `unikmer diff`: fix bug of hanging when the first file having no Kmers.

### [v0.4.2](https://github.com/shenwei356/unikmer/releases/tag/v0.4.2) - 2018-10-13 [![Github Releases (by Release)](https://img.shields.io/github/downloads/shenwei356/unikmer/v0.4.2/total.svg)](https://github.com/shenwei356/unikmer/releases/tag/v0.4.2)

- `unikmer stats/diff`: more intuitional output

### [v0.4.1](https://github.com/shenwei356/unikmer/releases/tag/v0.4.1) - 2018-10-10 [![Github Releases (by Release)](https://img.shields.io/github/downloads/shenwei356/unikmer/v0.4.1/total.svg)](https://github.com/shenwei356/unikmer/releases/tag/v0.4.1)

- Better performance of writing and reading binary files 

### [v0.4.0](https://github.com/shenwei356/unikmer/releases/tag/v0.4.0) - 2018-10-09 [![Github Releases (by Release)](https://img.shields.io/github/downloads/shenwei356/unikmer/v0.4.0/total.svg)](https://github.com/shenwei356/unikmer/releases/tag/v0.4.0)

- **Binary serialization format changed.**
- new command `unikmer sort`: sort binary files
- `unikmer count/diff/union/inter`: better performance, add option to sort Kmers which significantly reduces file size
- `unikmer dump`: changed option
- `unikmer count`: changed option

### [v0.3.1](https://github.com/shenwei356/unikmer/releases/tag/v0.3.1) - 2018-09-25 [![Github Releases (by Release)](https://img.shields.io/github/downloads/shenwei356/unikmer/v0.3.1/total.svg)](https://github.com/shenwei356/unikmer/releases/tag/v0.3.1)

- **Binary serialization format changed.**
- new command `unikmer stats`: statistics of binary files.
- `unikmer`: adding global option `-i/--infile-list` for reading files listed in file.
- `unikmer diff`: fixed a concurrency bug when no diff found.

### [v0.2.1](https://github.com/shenwei356/unikmer/releases/tag/v0.2.1) - 2018-09-23 [![Github Releases (by Release)](https://img.shields.io/github/downloads/shenwei356/unikmer/v0.2.1/total.svg)](https://github.com/shenwei356/unikmer/releases/tag/v0.2.1)

- `unikmer count`: performance improvement and new option `--canonical` for only keeping canonical Kmers.

### [v0.2.0](https://github.com/shenwei356/unikmer/releases/tag/v0.2) - 2018-09-09 [![Github Releases (by Release)](https://img.shields.io/github/downloads/shenwei356/unikmer/v0.2/total.svg)](https://github.com/shenwei356/unikmer/releases/tag/v0.2)

- new command `unikmer sample`: sample Kmers from binary files.
- new global options:
- `-c, --compact`:     write more compact binary file with little loss of speed.
- `-C, --no-compress`:   do not compress binary file (not recommended).
- some improvements.

### [v0.1.0](https://github.com/shenwei356/unikmer/releases/tag/v0.1.0) - 2018-08-09 [![Github Releases (by Release)](https://img.shields.io/github/downloads/shenwei356/unikmer/v0.1.0/total.svg)](https://github.com/shenwei356/unikmer/releases/tag/v0.1.0)

- first release
