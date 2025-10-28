#Set the parent directory (this should be the same directory as is set in the common_variables.sh script

# This file has the basic info used by multiple scripts in the subfolders handling the RPF and RNA reads/samples
parent_dir <- 'this/is/your/main/data/folder' #This is the path to the parent directory that contains all the data and where all the processed data will be saved - same as used for Shell scripts
# setwd(paste0(parent_dir, "/R_scripts")) # if needed

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

