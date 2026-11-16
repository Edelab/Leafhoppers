#!/bin/bash

ulimit -S -n 12288
Fspades="03_spades"
Fblast="04_blast"
Fmito_fragment="06_mitochondira_fragment"
Ffull_mito="07_full_mitochondria"

module load ncbiblast/2.16.0
for i in $(ls $Fspades)
do 
	grep -a mitoc $Fblast/blast_"$i"_contig | cut -f1 | uniq >"$Fblast"/"$i"_mitoseq.tmp
	samtools faidx -r "$Fblast"/"$i"_mitoseq.tmp $Fspades/$i/contigs.fasta >"$Fmito_fragments"/$i.fasta
	samtools faidx "$Fmito_fragments"/"$i".fasta
	awk '$2>15000 && $2<20000' "$Fmito_fragments"/$i.fasta.fai | cut -f 1 >complete_mito.tmp
	if [ -s complete_mito.tmp ]
	then
		sed -i "s/^/\^/" complete_mito.tmp
		sed -i "s/$/\t/" complete_mito.tmp
		grep -f complete_mito.tmp  $Fblast/blast_"$i"_contig | grep "Gonatopus" | cut -f 1 | sort | uniq >complete_gonatopus.tmp
		samtools faidx -r complete_gonatopus.tmp $Fspades/$i/contigs.fasta > $"Ffull_mito"/"$i".fasta
	fi

rm "$BLAST"/"$sample"_mitoseq.tmp complete_mito.tmp complete_gonatopus.tmp
done

