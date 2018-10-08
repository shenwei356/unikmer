#!/bin/sh
f=$1;

echo -e "K\tplain\tplain.gz\t.unik\t.unik.cpt\t.unik.sort\t.unik.ungz\t.unik.cpt.ungz\tplain.gz(%)\t.unik(%)\t.unik.cpt(%)\t.unik.sort(%)\t.unik.ungz(%)\t.unik.cpt.ungz(%)"
seq 1 2 32 \
    | rush -v f=$f -k ' \
        unikmer count -K -k {} {f} -o {f}.k{}; \
        unikmer count -K -k {} {f} -o {f}.k{}.c -c; \
        unikmer count -K -k {} {f} -o {f}.k{}.s -c -s; \
        unikmer count -K -k {} {f} -o {f}.k{}.C -C; \
        unikmer count -K -k {} {f} -o {f}.k{}.cC -c -C; \
        unikmer view {f}.k{}.unik > {f}.k{}.unik.plain; \
        gzip -c -k {f}.k{}.unik.plain > {f}.k{}.unik.plain.gz; \
          \
        p=$(ls -l {f}.k{}.unik.plain | cut -d " " -f 5); \
        gp=$(ls -l {f}.k{}.unik.plain.gz | cut -d " " -f 5); \
        o=$(ls -l {f}.k{}.unik | cut -d " " -f 5); \
        c=$(ls -l {f}.k{}.c.unik | cut -d " " -f 5); \
        s=$(ls -l {f}.k{}.s.unik | cut -d " " -f 5); \
        C=$(ls -l {f}.k{}.C.unik | cut -d " " -f 5); \
        cC=$(ls -l {f}.k{}.cC.unik | cut -d " " -f 5); \
        echo -e "{}\t$p\t$gp\t$o\t$c\t$s\t$C\t$cC" \
        ' \
    | csvtk -H -t mutate2 -L 1 -e '$3/$2*100' \
    | csvtk -H -t mutate2 -L 1 -e '$4/$2*100' \
    | csvtk -H -t mutate2 -L 1 -e '$5/$2*100' \
    | csvtk -H -t mutate2 -L 1 -e '$6/$2*100' \
    | csvtk -H -t mutate2 -L 1 -e '$7/$2*100' \
    | csvtk -H -t mutate2 -L 1 -e '$8/$2*100' \

ls $f.k* | rush '/bin/rm {}'
