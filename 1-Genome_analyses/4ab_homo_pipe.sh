#!/bin/bash

export PATH="$PATH:/usr_storage/jcf/.conda/envs/BUSCO"
source /usr_storage/jcf/geta-user204/.bashrc


rna_1_fq="cat $1|sed ":a;N;s/\n/,/g;ta"" #1.fq path list 
rna_2_fq="cat $2|sed ":a;N;s/\n/,/g;ta"" #2.fq path list
genome="$3" #genome fasta file 
conf="$4" #small genome conf.txt of geta pipepline setting as default parameters
out=${genome%.*}
homo_pro="$5"

geta.pl \
	--RM_species Embryophyta\
	--out_prefix `pwd`/$out \
	--config $conf \
	--cpu 20 \
	--protein $homo_pro\
	-genome $genome \
	-1 $rna_1_fq \
	-2 $rna_2_fq \
	--augustus_species $out

