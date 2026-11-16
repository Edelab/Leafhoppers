#!/bin/bash
# 16 CPUs
# 150 Go

# assemble genome from all trimmed dataset in 02_fastp


FASTP="02_fastp"
SPADES="03_spades"


module load SPAdes/4.1.0
# Iterate over files in data folder
for file in $(ls -1 $rawdir/*_1*.gz)
do

    input_file=$(echo "$file" | perl -pe 's/_1.*\.fastq.gz//')
    output_file=$(basename "$input_file")
    echo "Assembling: $output_file"
    spades.py -o $SPADES/$input_file \
	      -1 ${input_file}_1.fastq.gz \
	      -2 ${input_file}_2.fastq.gz \
	      --careful \
	      -t 16 \
	      -k 21,33,55,77,99,127
done

