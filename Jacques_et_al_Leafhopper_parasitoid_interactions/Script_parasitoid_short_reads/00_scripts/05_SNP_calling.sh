#!/bin/bash


module load bwa/0.7.17
module load SAMtools/1.18
module load java/jdk/1.8.0_102
module load bamUtil
module load bcftools/1.12

Reference_entry="PV999721.1"
curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=${Reference_entry=}&rettype=fasta&retmode=txt">../$Reference_entry.fasta


Freads="02_fastp"
Reference="../$Reference_entry.fasta"
Fout="03_bams"
Fvcf="06_vcf"
NCPU=10
MARKDUPS="/prg/picard-tools/1.119/MarkDuplicates.jar"
BAMINDEX="/prg/picard-tools/1.119/BuildBamIndex.jar"
CSD="/prg/picard-tools/1.119/CreateSequenceDictionary.jar"
GATK="/home/clrou103/00-soft/GATK/GenomeAnalysisTK.jar"


bwa index $Reference
samtools faidx $Reference
java -jar $CDS R=$Reference O="${Reference_entry/.fasta}.dict"


# Iterate over files in data folder
for file in $(ls -1 $Freads/*_1*.gz)
do
	input_file=$(echo "$file" | perl -pe 's/_1.*\.fastq.gz//')
    	output_file=$(basename "$input_file")
        ID="@RG\tID:$input_file\tSM:$input_file\tPL:Illumina"
	
	#align reads
	bwa mem -t $NCPU -R $ID $Reference "$input_file"_1.fastq.gz "$input_file"_2.fastq.gz | 
	samtools view -Sb - >$Fout/$output_file.bam
	
	# Sort reads
	samtools sorts -@ $NCPU $Fout/$output_file.bam > $Fout/$output_file.sorted.bam
	samtools index $Fout/$output_file.sorted.bam
	rm $Fout/$output_file.bam

	# Mark Duplicates
	java -jar "$MARKDUPS" INPUT=$Fout/$output_file.sorted.bam \
	       		      OUTPUT="$Fout"/$output_file.dedup.bam \
			      METRICS_FILE=$Fout/$output_file.metrics.txt  \
			      VALIDATION_STRINGENCY=SILENT \
			      REMOVE_DUPLICATES=true \
			      MAX_FILE_HANDLES=100 \
			      TMP_DIR="./tmp"

	rm $Fout/$output_file.sorted.bam

	# index BAM
	java -jar $BAMINDEX INPUT="$Fout"/$output_file.dedup.bam

	# Realign Indels
	java -jar "$GATK" -T RealignerTargetCreator -R $Reference -I "$Fout"/$output_file.dedup.bam -o "$Fout"/$output_file.dedup.bam.intervals
	java -jar $GATK -T IndelRealigner -R $Reference \
		-I "$Fout"/$output_file.dedup.bam \
		-targetIntervals "$Fout"/$output_file.dedup.bam.intervals \
		--consensusDeterminationModel USE_READS  \
		-o "$Fout"/$output_file.realigned.bam
	
	rm "$Fout"/$output_file.dedup.bam "$Fout"/$output_file.dedup.bam.intervals
	
	#ClipOverlap
	bam clipOverlap --in "$Fout"/$output_file.realigned.bam \
        		--out "$Fout"/$output_file.temp.bam \
        		--unmapped --storeOrig OC --stats
	samtools view -h -F 4 "$Fout"/$output_file.temp.bam |\
		samtools view -hb -e "sclen <30" \ > "$Fout"/$output_file.final.bam
	samtools index "$Fout"/$output_file.final.bam
	rm "$Fout"/$output_file.temp.bam


done



# SNP calling
ls $Fout/*.final.bam > list_bams
bcftools mpileup -f $Reference -a AD,DP,SP,ADF,ADR -q 5 -Ou -d 20000 $list_bam |\
	bcftools call -m -a GP,GQ -Oz -G - -o $Fvcf/raw.vcf.gz --ploidy 1

python 00_utils/Filtering_coverage_ratio.py --vcf $Fvcf/raw.vcf.gz -lratio 0.1 -hratio 0.9

