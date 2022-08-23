#!/bin/bash


export PATH="$PATH:/usr_storage/jcf/.conda/envs/PASA "
source  /pub_storage2/new_PASA/.bashrc

genome="$1" #genome fasta file
annotation_conf="$2" #pasa annotation compare conf 
transcripts_fasta="$3" #transcripts_fasta file for PASA seqclean step
gff3="$4" #gff3 for PASA updata


Launch_PASA_pipeline.pl \
		-c $annotation_conf\
		-A -T -L \
		-g $genome\
		-t ${transcripts_fasta}.clean \
		-u $transcripts_fasta \
		--annots $gff3