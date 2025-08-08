#!/usr/bin/env bash

# This scripts extracts read counts from the log files generated at the processing steps in order to calculate later on the final read depths per sample
# These are used in QC steps in R

#read in variables
source common_variables.sh

#monosomes
for filename in $Totals_filenames
do
extract_read_counts.py ${filename} Totals -log_dir $log_dir &
done
wait
