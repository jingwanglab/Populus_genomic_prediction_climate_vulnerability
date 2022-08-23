#!/bin/bash

export PATH="$PATH:/usr_storage/jcf/.conda/envs/PASA"
source  /pub_storage2/new_PASA/.bashrc

#cat $Trinity_GG $Trinity_denovo >transcripts.fasta #
transcripts_fasta="$1" # transcripts.fasta generated from merging fasta file of Trinity denovo and Trinity genome guided mode

#perl -e 'while(<>) { print "$1\n" if />(\S+)/ }' Trinity.fasta >tdn.accs #
denovo_transcript_id="$2" 
alignAssembly_config="$3"
genome="$4" #reference fasta file



seqclean $transcripts_fasta \
	     -v /pub_storage2/PASA/UniVec
		 

Launch_PASA_pipeline.pl -c $alignAssembly_config \
					    -C -R -T \
						-g $genome \ 
						-t $transcripts_fasta.clean \
						-u ${transcripts_fasta} \
						--ALIGNERS gmap,blat \
						--CPU 8 \ 
						--TDN $denovo_transcript_id
