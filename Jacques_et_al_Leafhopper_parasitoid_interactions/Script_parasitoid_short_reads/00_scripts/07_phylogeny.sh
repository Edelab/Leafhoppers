#!/bin/bash


# Fetching accession fasta
for i in PP738850 PP738821 PP738805 PP600749 MF905225 MZ628882 MZ628779 MZ628513 MZ628452 MZ627964 MZ627864 MZ627284 MZ625137 MW983030 MW981730 JN306858 JN306857 MG377583;
do 
	curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=${i}&rettype=fasta&retmode=txt">>08_NCBI_fasta/accession.fasta;
       	sleep 20;
done


#Cleaning names and empty lines
awk 'NF' 08_NCBI_fasta/accession.fasta >08_NCBI_fasta/accession_fasta.tmp
sed "s/ /_/" 08_NCBI_fasta/accession_fasta.tmp  | sed "s/>.*_/>/" | sed "s/ /_/"| sed "s/ .*//" > 08_NCBI_fasta/accession.fasta



# Alignement
cat 07_fastas/COI.fasta 08_NCBI_fasta/accession.fasta >09_alignement/COI.fasta
seqkit grep -v LR02_0002 09_alignement/COI.fasta > 09_alignement/Clean_COI.fasta
# Alignment step was performed in genenious using default alignement parameters and Geneious algorithm
# You need to import 09_alignement/Clean_COI.fasta then export the resulting alignemnt file as fasta in 09_alignement/Clean_COI.ali

# We then used Gblock with default parameters to clean using default value, which creates the file
# 09_alignement/Clean_COI.ali.gb.fa which we used to build the tree


# IqTree
iqtree -s 09_alignement/Clean_COI.ali.gb.fa -st DNA -m MF -B 1000




# IqTree for Multiple genes:
seqkit concat -o 09_alignement/ALL_gene.ali 07_fastas/cox1.fasta 07_fastas/atp6.fasta 07_fastas/cob.fasta 07_fastas/cox2.fasta 07_fastas/cox3.fasta 07_fastas/nad4.fasta
seqkit grep -v LR02_0002 09_alignement/ALL_gene.ali 09_alignement/Clean_ALL_gene.ali
# We then used Gblock with default parameters to clean using default value, which creates the file
# 09_alignement/Clean_ALL_gene.ali.gb.fa which we used to build the tree

#IqTree
iqtree -s 09_alignement/Clean_ALL_gene.ali.gb.fa -st DNA -m MFP -B 1000

