#!/bin/bash


module load SAMtools/1.18
module load bcftools/1.12

Reference_entry="PV999721.1"
curl "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=gff3&id=PV999721.1" >../$Reference_entry.gff


Freads="02_fastp"
Reference="../$Reference_entry.fasta"
GFF="../$Reference_entry.fasta"
Ffasta="07_fastas"
Fvcf="06_vcf"


# Buildign consensus
cat $GFF| \
grep -P "\tgene\t"| \
while read line
    do
        gene=$(echo $line |cut -d" " -f9  |sed "s/;.*Name=//" |sed "s/;.*//")
        pos=$(echo $line | cut -d" " -f1,4,5 | sed "s/ /:/"|sed "s/ /-/")

        for file in $(ls -1 $Freads/*_1*.gz)
        do
                input=$(echo "$file" | perl -pe 's/_1.*\.fastq.gz//')
                SN=$(basename $input)

                echo ">$file" >>$Ffasta/"$gene".fasta
                samtools faidx $Reference $pos | bcftools consensus -a N -M N -s $SN $Fvcf/ratio_qual_filtered_0.1_0.9.vcf | \
                grep -v ">" >>$Ffasta/"$gene".fasta
            done

    done

