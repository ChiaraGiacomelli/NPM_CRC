#!/usr/bin/env bash

# variables and folders to be used for Ribominer package
# put the genome_sorted.bam files also in the RiboMiner folder
# from https://github.com/xryanglab/RiboMiner

parent_dir='this/is/your/main/data/folder' #This is the path to the parent directory that contains all the data and where all the processed data will be saved
fasta_dir='this/is/the/folder/with/all/your/organism/FASTAs'

STAR_GTF=${fasta_dir}/GENCODE/vM27/original/gencode.vM27.annotation.gtf
PRIMARY_ASSEMBLY=${fasta_dir}/GENCODE/vM27/original/GRCm39.primary_assembly.genome.fa

RM_ATTRIBUTES=${parent_dir}/RiboMiner/attributes.txt

cd ../RiboMiner

# Gets coordinates for all protein coding genes & the fasta with all the sequences - puts them in the RiboMiner folder
prepare_transcripts -g $STAR_GTF -f $PRIMARY_ASSEMBLY -o $parent_dir/RiboMiner

# Prepares the annotation file for the longest transcript per gene
OutputTranscriptInfo -c transcripts_cds.txt -g $STAR_GTF -f transcripts_sequence.fa -o longest.transcripts.info.txt -O all.transcripts.info.txt

# Prepares the sequence file for the longest transcript per gene
GetProteinCodingSequence -i transcripts_sequence.fa  -c longest.transcripts.info.txt -o NPM_KO --mode whole

# GitHub also suggests to get the UTR sequences
GetUTRSequences -i transcripts_sequence.fa -o NPM_KO -c transcripts_cds.txt

# Polarity score calculation
PolarityCalculation -f attributes.txt -c longest.transcripts.info.txt -o NPM_KO -n 64
