#!/usr/bin/env bash

# At the end of this script, the files are to be explored with R scripts to plot for features such as periodicity, lengths, localisation on the CDS vs UTRs etc
# Once the Plots have been generated and checked, proceed to the Step F where the sizes and correct offsets are used for counting the final reads

#read in variables
source common_variables.sh

#make an fai (fasta index) file from the fasta using samtools. This is required for the counting script and needs to exist before running counting_script.py
samtools faidx $most_abundant_fasta

#run the counting_script.py with a range of read lengths (adjust below if required, currently set to 25-35)
for filename in $RPF_filenames
do
for length in $(seq 27 38)
do
counting_script.py -bam $BAM_dir/${filename}_pc_deduplicated_sorted.bam -fasta $most_abundant_fasta -len $length -out_file ${filename}_pc_L${length}_Off0.counts -out_dir $counts_dir &
done
done
wait

#set offset
offset=15

#run summing_region_counts.py script
for filename in $RPF_filenames
do
for length in $(seq 27 38)
do
summing_region_counts.py ${filename}_pc_L${length}_Off0.counts $offset $region_lengths -in_dir $counts_dir -out_dir $region_counts_dir &
done
done
wait

#set number of nt to splice
n=50

#run summing_spliced_counts.py script
for filename in $RPF_filenames
do
for length in $(seq 27 38)
do
summing_spliced_counts.py ${filename}_pc_L${length}_Off0.counts $n $region_lengths -in_dir $counts_dir -out_dir $spliced_counts_dir &
done
done
wait

#run periodicity.py script
for filename in $RPF_filenames
do
for length in $(seq 27 38)
do
periodicity.py ${filename}_pc_L${length}_Off0.counts $region_lengths -offset $offset -in_dir $counts_dir -out_dir $periodicity_dir &
done
done
wait
