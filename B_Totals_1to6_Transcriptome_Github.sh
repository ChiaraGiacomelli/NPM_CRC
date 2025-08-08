#!/usr/bin/env bash

#This script runs cutadapt on all your raw fastq files. For more info on settings see https://cutadapt.readthedocs.io/en/stable/guide.html

#You need to make sure the 3' adaptor sequence in the variables.txt file is correct. The -a option specifies this here

#The --nextseq-trim=20 option will trim bases from the 3' end if the quality score is below 20. This is only for sequencing data generated on the
#Next-seq (which is what we have here at the Beatson). If using external data set that has been sequenced on a Hi-seq, replace this command with -q 20
#It should state which sequencing platform was used on the GEO page
#The following text from the cutadpat manual explains why this is

#Some Illumina instruments use a two-color chemistry to encode the four bases. This includes the NextSeq and the NovaSeq.
#In those instruments, a ‘dark cycle’ (with no detected color) encodes a G. However, dark cycles also occur when sequencing “falls off” the end of the fragment.
#The read then contains a run of high-quality, but incorrect “G” calls at its 3’ end.
#Since the regular quality-trimming algorithm cannot deal with this situation, you need to use the --nextseq-trim option:
#This works like regular quality trimming (where one would use -q 20 instead), except that the qualities of G bases are ignored.

#-m specifies the minimum read lengths following adaptor removal and base trimming, which we set as 30. As the fragment length for total RNA-seq is normally longer than the sequencing length no upper limit is specified

#1> causes all the text that is normally printed to the screen to be saved in a log file in your logs directory for each sample

#once cutadapt is complete, fastQC is ran on the output fastq files to check they are as expected

#read in variables
source common_variables.sh

# CUTADAPTOR
#run cutadapt
for filename in $Totals_filenames
do
cutadapt $fastq_dir/${filename}.fastq -a $Totals_adaptor --nextseq-trim=20 -m 30 --cores=0 -o $fastq_dir/${filename}_cutadapt.fastq 1> $log_dir/${filename}_cutadapt_log.txt
done

#run fastqc on cutadapt output
for filename in $Totals_filenames
do
fastqc $fastq_dir/${filename}_cutadapt.fastq --outdir=$fastqc_dir &
done
wait

# DEDUPLICATE STEP 1 - umi removal
#read deduplication
for filename in $Totals_filenames
do
umi_tools extract -I $fastq_dir/${filename}_cutadapt.fastq -S $fastq_dir/${filename}_UMI_clipped.fastq --bc-pattern=NNNNNNNNNNNN --log=$log_dir/${filename}_extracted_UMIs.log &
done
wait

#run fastqc on output
for filename in $Totals_filenames
do
fastqc $fastq_dir/${filename}_UMI_clipped.fastq --outdir=$fastqc_dir &
done
wait

# ALIGN
#Align to protein coding transcriptome
for filename in $Totals_filenames
do
bowtie2 -S $SAM_dir/${filename}_pc.sam -U $fastq_dir/${filename}_UMI_clipped.fastq -x $rsem_index --threads $threadN --sensitive --dpad 0 --gbar 99999999 --mp 1,1 --np 1 --score-min L,0,-0.1 2> $log_dir/${filename}_pc_log.txt
done

#convert sam to bam 
for filename in $Totals_filenames
do
samtools view -b $SAM_dir/${filename}_pc.sam > $BAM_dir/${filename}_pc.bam &
done
wait

#sort bam
for filename in $Totals_filenames
do
samtools sort $BAM_dir/${filename}_pc.bam -o $BAM_dir/${filename}_pc_sorted.bam -@ $threadN -m 1G
done

#index bam
for filename in $Totals_filenames
do
samtools index $BAM_dir/${filename}_pc_sorted.bam $BAM_dir/${filename}_pc_sorted.bai &
done
wait

# DEDUPLICATE STEP 2
#run UMI tools deduplication function
for filename in $Totals_filenames
do
umi_tools dedup -I $BAM_dir/${filename}_pc_sorted.bam -S $BAM_dir/${filename}_pc_deduplicated.bam --output-stats=$log_dir/${filename}_deduplication 1> $log_dir/${filename}_deduplication_log.txt &
done
wait

# Isoform quantification
#Align to protein coding transcriptome
for filename in $Totals_filenames
do
rsem-calculate-expression --strandedness forward --fragment-length-mean 300 --fragment-length-sd 100 --alignments $BAM_dir/${filename}_pc_deduplicated.bam $rsem_index $rsem_dir/${filename} &
done
wait

# Write MAT
filter_FASTA.py $pc_fasta $most_abundant_transcripts_dir/most_abundant_transcripts.txt $most_abundant_fasta
