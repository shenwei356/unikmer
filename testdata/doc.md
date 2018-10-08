## run

    f=Ecoli-MG1655.fasta.gz
    # f=A.muciniphila-ATCC_BAA-835.fasta.gz
    # f=Leishmania-donovani_BPK282A1.fasta.gz # FR799588.2
    ./cr.sh $f > $f.cr.tsv

    # plot
    csvtk cut -t -f 1,9-14 $f.cr.tsv \
        | csvtk -t rename2 -f 2-7 -p "\(%\)" \
        | csvtk -t gather  -k group  -v value -f 2-7 \
        | csvtk plot line -t -x 1 -y 3 -g 2 \
            --x-min 3 --x-max 45  --y-max 100 --ylab "compression rate (%)" \
            --width 4.5 --height 3.3 \
            --title $f \
        > $f.cr.tsv.png 
