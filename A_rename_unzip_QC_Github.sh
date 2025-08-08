#!/usr/bin/env bash

#read in variables
source common_variables.sh

# remove name parts from file name
for f in $fastq_dir/*.fastq.gz; do newname=$( echo $f | sed -r 's/_S[0-9]{2}_R1_001//' ); mv $f $newname; done
for f in $fastq_dir/*.fastq.gz; do newname=$( echo $f | sed -r 's/_S[0-9]{1}_R1_001//' ); mv $f $newname; done

#unzip
for filename in $RPF_filenames
do
gunzip $fastq_dir/${filename}.fastq.gz &
done
wait

for filename in $Totals_filenames
do
gunzip $fastq_dir/${filename}.fastq.gz &
done
wait

#run fastQC on RPFs
for filename in $RPF_filenames
do
fastqc $fastq_dir/${filename}.fastq --outdir=$fastqc_dir &
done
wait

#run fastQC on Totals
for filename in $Totals_filenames
do
fastqc $fastq_dir/${filename}.fastq --outdir=$fastqc_dir &
done
wait
