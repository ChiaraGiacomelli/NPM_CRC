#!/usr/bin/env bash

# This scripts extracts read counts from the log files generated at the processing steps in order to calculate later on the final read depths per sample
# These are used in QC steps in R

#read in variables
source common_variables.sh

for filename in $RPF_filenames
do
extract_read_counts.py ${filename} RPFs -log_dir $log_dir &
done
wait
