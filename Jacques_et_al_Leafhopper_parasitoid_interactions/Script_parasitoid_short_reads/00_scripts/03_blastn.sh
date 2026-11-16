#!/bin/bash

ulimit -S -n 12288
Fspades="03_spades"
Fblast="04_blast"


module load ncbiblast/2.16.0
for i in $(ls $Fspades)
do 
	srun -c 10 --mem 150G -p medium --time 07-00:00:00 -o blast_$i.log \
	       	blastn -db /biodata/blast/nt \
		-query $Fspades/$i/contigs.fasta \
		-out $Fblast/blast_"$i"_contig \
		-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle" \
	       	-max_target_seqs 10 -num_threads 10& 

done

