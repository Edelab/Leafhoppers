#!/bin/bash
# 6 CPUs
# 10 Go

# trim input files in 04_raw_data with fastp
# No need to parallelize more, using 6 CPUs makes fastp really fast
# Adapted from https://github.com/enormandeau/wgs_sample_preparation
# WARNING: File name formats may require changing code below
rawdir=$1





echo -e ">TrueSeq R1\nAGATCGGAAGAGCACACGTCTGAACTCCAGTCA\n>TrueSeq R2\nAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT\n>Reverse Trueseq R1\nTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT\n>Reverse Trueseq R2\nACACTCTTTCCCTACACGACGCTCTTCCGATCT">adatapters.fa

module load fastp
# Iterate over files in data folder
for file in $(ls -1 $rawdir/*_R1*.gz)
do
    input_file=$(echo "$file" | perl -pe 's/_R1.*\.fastq.gz//')
    output_file=$(basename "$input_file")
    echo "Treating: $output_file"

    fastp -w 6 \
        -i "$input_file"_R1.fastq.gz \
        -I "$input_file"_R2.fastq.gz \
        -o 02_fastp/"$output_file"_1.trimmed.fastq.gz \
        -O 02_fastp/"$output_file"_2.trimmed.fastq.gz \
        -j 02_faspt/01_reports/"$output_file".json \
        -h 02_fastp/01_reports/"$output_file".html \
        --adapter_fasta adatapters.fa
done

