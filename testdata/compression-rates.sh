#!/bin/sh
f=$1;

echo -e "K\tplain\t.gz\t.unik\t.gz(%)\t.unik(%)"
seq 1 2 31 \
    | rush -v f=$f -k \
        'unikmer count -k {} {f} -o {f}.k{}; \
        unikmer view {f}.k{}.unik > {f}.k{}.unik.plain; \
        gzip -c -k {f}.k{}.unik.plain > {f}.k{}.unik.plain.gz; \
        a=$(ls -l {f}.k{}.unik.plain | cut -d " " -f 5); \
        b=$(ls -l {f}.k{}.unik | cut -d " " -f 5); \
        g=$(ls -l {f}.k{}.unik.plain.gz | cut -d " " -f 5); \
        echo -e "{}\t$a\t$g\t$b"' \
    | csvtk -H -t mutate2 -L 1 -e '$3/$2*100' \
    | csvtk -H -t mutate2 -L 1 -e '$4/$2*100'

ls $f.k* | rush 'rm {}'
