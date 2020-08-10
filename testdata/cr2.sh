#!/bin/sh
f=$1;

function NK (){
    echo $(($1<<10))
}

function NM (){
    echo $(($1<<20))
}

function NG (){
    echo $(($1<<30))
}


# echo counting k-mers (not sorted)
j=12 # depends on your cpus and memory 
# time seq 11 2 32  \
#     | rush -j $j -v f=$f 'if [[ -e {f}.k{}.unik ]]; then exit; fi; \
#         unikmer count -K -k {} {f} -o {f}.k{}'  

 
echo -e  "k\tnum\tplain\tgzip\tunik\tcunik\tsunik\tunikC\tcunikC\tsunikC";
time for s in $(NK 1) $(NM 1) $(NM 10) $(NM 100); do
    seq 11 2 32 \
        | rush -j $j -k -v f=$f -v s=$s -v 'F={f}.k{}.unik' \
            'sum=$(unikmer stats {F} -Ta | csvtk cut -t -f number | sed 1d);
            if [[ $sum -lt {s} ]]; then exit; fi; \
            plain=$(unikmer head -n {s} {F} -C | unikmer view | wc -c); \
             gzip=$(unikmer head -n {s} {F} -C | unikmer view | gzip -c | wc -c); \
             unik=$(unikmer head -n {s} {F} | wc -c); \
            cunik=$(unikmer head -n {s} {F} | unikmer view | unikmer dump -c | wc -c); \
            sunik=$(unikmer head -n {s} {F} | unikmer sort | wc -c); \
            unikC=$(unikmer head -n {s} {F} -C | wc -c); \
           cunikC=$(unikmer head -n {s} {F} -C | unikmer view -C | unikmer dump -c -C | wc -c); \
           sunikC=$(unikmer head -n {s} {F} -C | unikmer sort -C | wc -c); \
            echo -ne "{}\t{s}\t$plain\t$gzip\t$unik\t$cunik\t$sunik\t$unikC\t$cunikC\t$sunikC\n";'
done

