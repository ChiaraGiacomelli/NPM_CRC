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

#-m and -M specify the minimum and Maximum read lengths following adaptor removal and base trimming. RPFs should be roughly 30nt but can vary.
#However if the libraries also contain UMIs then you need to add this on too. UMIs from the nextflex kit that we use are 4nt on each end of the read
#therefore a 30nt RPF should be 38nt with UMIs following adaptor removal. Using -m 30 -M 50 means after UMI removal you will be left with read lengths 22-42
#which is suitable for this situation but will need to be modified to suit specific needs if using external data which may have different UMI lengths or may not contain UMIs

#--cores=0 will automatically detect the number of cores available on the system and use this amount

#1> causes all the text that is normally printed to the screen to be saved in a log file in your logs directory for each sample

#once cutadapt is complete, fastQC is ran on the output fastq files to check they are as expected

#read in variables
source common_variables.sh

#run cutadapt
for filename in $RPF_filenames
do
cutadapt $fastq_dir/${filename}.fastq -a $RPF_adaptor --nextseq-trim=20 -m 30 -M 50 --cores=0 -o $fastq_dir/${filename}_cutadapt.fastq 1> $log_dir/${filename}_cutadapt_log.txt
done

#run fastqc on cutadapt output
for filename in $RPF_filenames
do
fastqc $fastq_dir/${filename}_cutadapt.fastq --outdir=$fastqc_dir &
done
wait

# DEDUPLICATION STEP 1
#extract UMIs
for filename in $RPF_filenames
do
umi_tools extract -I $fastq_dir/${filename}_cutadapt.fastq --extract-method=regex --bc-pattern='^(?P<umi_1>.{4}).+(?P<umi_2>.{4})$' -S $fastq_dir/${filename}_UMI_clipped.fastq --log=$log_dir/${filename}_extracted_UMIs.log &
done
wait

#run fastqc on output
for filename in $RPF_filenames
do
fastqc $fastq_dir/${filename}_UMI_clipped.fastq --outdir=$fastqc_dir &
done
wait

# Alignments
#Align to rRNA
for filename in $RPF_filenames
do
bbmap.sh in=$fastq_dir/${filename}_UMI_clipped.fastq ref=$rRNA_fasta outm=$fastq_dir/${filename}_rRNA.fastq outu=$fastq_dir/${filename}_non_rRNA.fastq ambiguous=best nodisk threads=$threadN 2> $log_dir/${filename}_rRNA_log.txt
done

#Align to tRNA fasta
for filename in $RPF_filenames
do
bbmap.sh in=$fastq_dir/${filename}_non_rRNA.fastq ref=$tRNA_fasta outm=$fastq_dir/${filename}_tRNA.fastq outu=$fastq_dir/${filename}_non_rRNA_tRNA.fastq ambiguous=best nodisk threads=$threadN 2> $log_dir/${filename}_tRNA_log.txt
done

#Align to protein coding transcriptome
for filename in $RPF_filenames
do
bbmap.sh in=$fastq_dir/${filename}_non_rRNA_tRNA.fastq out=$BAM_dir/${filename}_pc.bam ref=$most_abundant_fasta outm=$fastq_dir/${filename}_pc.fastq outu=$fastq_dir/${filename}_unaligned.fastq ambiguous=best nodisk trimreaddescription=t threads=$threadN 2> $log_dir/${filename}_pc_log.txt
done

#sort bam
for filename in $RPF_filenames
do
samtools sort $BAM_dir/${filename}_pc.bam -o $BAM_dir/${filename}_pc_sorted.bam -@ $threadN -m 1G
done

#index bam
for filename in $RPF_filenames
do
samtools index $BAM_dir/${filename}_pc_sorted.bam $BAM_dir/${filename}_pc_sorted.bai &
done
wait

#run fastqc on mapped reads
for filename in $RPF_filenames
do
fastqc $fastq_dir/${filename}_rRNA.fastq --outdir=$fastqc_dir &
fastqc $fastq_dir/${filename}_tRNA.fastq --outdir=$fastqc_dir &
fastqc $fastq_dir/${filename}_pc.fastq --outdir=$fastqc_dir &
fastqc $fastq_dir/${filename}_unaligned.fastq --outdir=$fastqc_dir &
done
wait

# DEDUPLICATION STEP 2
#run UMI tools deduplication function
for filename in $RPF_filenames
do
umi_tools dedup -I $BAM_dir/${filename}_pc_sorted.bam -S $BAM_dir/${filename}_pc_deduplicated.bam --output-stats=$log_dir/${filename}_deduplication 1> $log_dir/${filename}_deduplication_log.txt &
done
wait

#sort bam
for filename in $RPF_filenames
do
samtools sort $BAM_dir/${filename}_pc_deduplicated.bam -o $BAM_dir/${filename}_pc_deduplicated_sorted.bam -@ $threadN -m 1G
done

#index bam
for filename in $RPF_filenames
do
samtools index $BAM_dir/${filename}_pc_deduplicated_sorted.bam $BAM_dir/${filename}_pc_deduplicated_sorted.bai &
done
wait

# At this stage, time to check for counts and sizes and then extracting counts of lengths required
