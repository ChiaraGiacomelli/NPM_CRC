#Set the parent directory (this should be the same directory as is set in the common_variables.sh script

# This file has the basic info used by multiple scripts in the subfolders handling the RPF and RNA reads/samples

# choose your Machine
# VM direct runs
#machine_dir <- '~/data'
# Rstudio on Windows at work
#machine_dir <- 'N:'
# VPN connected on Mac
machine_dir <- '/Volumes/data-1'

parent_dir <- paste0(machine_dir, '/CGIACOME/AAJ_NPM/20230622_RiboSeq_APCKRAS/Ribo-seq-Ribo-seq2.0')

# setwd(paste0(parent_dir, "/R_scripts"))

#set sample names
RPF_sample_names <- c('NPM_WT_1_RPFs', 'NPM_WT_2_RPFs', 'NPM_WT_3_RPFs', 'NPM_WT_4_RPFs', 'NPM_KO_1_RPFs', 'NPM_KO_2_RPFs', 'NPM_KO_3_RPFs', 'NPM_KO_4_RPFs')
Total_sample_names <- c('NPM_WT_1_Totals', 'NPM_WT_2_Totals', 'NPM_WT_3_Totals', 'NPM_WT_4_Totals', 'NPM_KO_1_Totals', 'NPM_KO_2_Totals', 'NPM_KO_3_Totals', 'NPM_KO_4_Totals')

RPF_sample_info <- data.frame(sample = RPF_sample_names,
                             condition = c(rep("WT", 4), rep("KO", 4)),
                             replicate = factor(rep(c("1", "2", "3", "4"), 2)))

Total_sample_info <- data.frame(sample = Total_sample_names,
                                condition = c(rep("WT", 4), rep("KO", 4)),
                                replicate = factor(rep(c("1", "2", "3", "4"), 2)))
NPM_cols <- c("WT" = "grey34", "KO" = "red3")
WT_color = "grey34"
KO_color = "red3"

