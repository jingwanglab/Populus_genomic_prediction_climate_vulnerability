#!/bin/bash

export PATH="/usr_storage/xyf/jcf/genewise/EVM/EVidenceModeler-1.1.1/EvmUtils/:$PATH"

genome="$1" #genome fasta file 
augustus_gff3="$2" #gff3 generated from augutus 
genewise_gff3="$3" #gff3 generated from tblastn and genewise
pasa_align_gff3="$4" #gff3 generated from PASA 
repeat_gff3="$5" #repeat gff3 generated from repeatemasker
partition="$6" #partition path for evm



partition_EVM_inputs.pl \
		--genome $genome\
		--gene_predictions $augustus_gff3 \
		--protein_alignments $genewise_gff3 \
		--transcript_alignments $pasa_align_gff3 \
		--repeats $repeat_gff3 \
		--segmentSize 5000000 \
		--overlapSize 10000 \
		--partition_listing $partition
		
write_EVM_commands.pl \
		--genome $genome \
		--gene_predictions $augustus_gff3 \
		--protein_alignments $genewise_gff3 \
		--transcript_alignments $pasa_align_gff3 \
		--repeats $repeat_gff3 \
		--output_file_name evm.out \
		--weights $weight >command.list
		
ParaFly -c command.list -CPU 32 

recombine_EVM_partial_outputs.pl \
		--partitions $partition \
		--output_file_name evm.out 
		
convert_EVM_outputs_to_GFF3.pl \
		--partitions $partition \
		--output_file_name evm.out \
		--genome  $genome 

cat */evm.out.gff3 >evm.out.gff3