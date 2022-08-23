#!/bin/bash

genome=$1
index=${genome%.*}
rna_1_fq=`cat $2|grep 1P|sed ":a;N;s/\n/,/g;ta"` #1.fq path list
rna_2_fq=`cat $2|grep 2P|sed ":a;N;s/\n/,/g;ta"` #2.fq path list

#echo $index
hisat2-build -p 20 $genome $index

hisat2 -x $index \
           -1 $rna_1_fq\
           -2 $rna_2_fq\
           --threads 20 \
           --min-intronlen 20 \
           --max-intronlen 20000 \
           --dta \
           --score-min L,0.0,-0.4 \
           -S ${index}.sam


samtools sort -@ 20 \
                  -o ${index}.sorted.bam \
                      -O BAM \
                ${index}.sam
