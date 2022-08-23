#!/bin/bash

#conda activate trinity

export PATH="$PATH:/usr_storage/jcf/.conda/envs/trinity"

rna_1_fq="cat $1|sed ":a;N;s/\n/,/g;ta"" #1.fq path list 
rna_2_fq="cat $2|sed ":a;N;s/\n/,/g;ta"" #2.fq path list
bam="$3"  #sorted.bam from hisat
out=${bam%.*}


Trinity --left $rna_1_fq \
	    --right $rna_2_fq \
		--seqType fq  \
		--max_memory 100G \
		--no_normalize_reads \
		--CPU 20 \
		--bflyCalculateCPU  \
		--output trinity_denovo_$out
		
Trinity --genome_guided_bam $bam  \
		--genome_guided_max_intron 10000 \
		--max_memory 100G \
		--no_normalize_reads \
		--CPU 20 \
		--bflyCalculateCPU\
		--output trinity_GG_$out