## Compression rate comparison

No Taxids stored.

1. Prepare a genome sequence, I used human genome chromosome X (`t_chrX.fa.gz`)

2. Computation

    f=t_chrX.fa.gz
    ./cr2.sh $f > table.tsv
    
3. Plot
    
    
    cat table.tsv \
        | csvtk -t mutate2 -L 1 -n r_gzip -e '$gzip/$plain*100' \
        | csvtk -t mutate2 -L 1 -n r_unik.default -e '$unik/$plain*100' \
        | csvtk -t mutate2 -L 1 -n r_unik.compact -e '$cunik/$plain*100' \
        | csvtk -t mutate2 -L 1 -n r_unik.sorted -e '$sunik/$plain*100' \
        | csvtk -t cut -F -f k,num,r_* \
        | csvtk -t gather -k group -v value -F -f 'r_*' \
        | csvtk -t replace -f group -p 'r_' \
        | csvtk -t replace -f num -p '^(.+)$' -k size.tsv -r '{kv} k-mers' \
        > table.r.tsv
    
    ./plot.R
    
