# This script was developed by Dr.Pauline Herviou, Postdoctoral Scientist in the lab of Prof. Martin Bushell
# At the CRUK Scotland Institute, Glasgow

# The aim is to identify pause sites from Riboseq results

#load packages----
library(gridExtra)
library(tidyverse)

#read in common variables
source("common_variables.R")

#read in data----
load(file = file.path(parent_dir, "counts_list.Rdata"))

#import region lengths file
region_lengths <- read_csv(file = "/FASTAs/mouse/GENCODE/vM27/transcript_info/gencode.vM27.pc_transcripts_region_lengths.csv", col_names = c("transcript", "UTR5_len", "CDS_len", "UTR3_len"))

#import most abundant transcripts info
most_abundant_transcripts <- read_csv(file = file.path(parent_dir, "most_abundant_transcripts_IDs.csv"))

#data annotation function
annotate_data <- function(df, sampl) {
  df %>%
    select(c("transcript", "Position", "Nucleotide", "normalised_CPM")) %>% #select needed columns
    mutate(sample = sampl) %>%
    inner_join(region_lengths, by = "transcript") %>%
    mutate(region = factor(case_when(Position <= UTR5_len ~ 'utr5', #annotate transcript regions
                                   Position <= UTR5_len + CDS_len ~ 'cds',
                                   TRUE ~ 'utr3')),
         cds_position = Position - UTR5_len) %>% #calculate CDS position
  filter(region == "cds") %>% #keep only CDS data
  select(-any_of(c("UTR5_len", "CDS_len", "UTR3_len", "region"))) -> annotated_df
  
  return(annotated_df)
}  

#Select CDS positions and select needed columns only for all conditions
annotated_list <- imap(counts_list, annotate_data) #compared to lapply, imap inputs object names/indexes (here sample names as "sampl") into the function
print(head(annotated_list[[1]]))

#concatenate tables from each sample 
all_counts <- do.call("rbind", annotated_list)

#get sample info data from sample names
all_counts %>%
  separate_wider_delim(sample, delim="_", names=c("npm", "condition", "replicate", "seq")) %>%
  select(-any_of(c("seq", "npm"))) -> all_mean_counts

print(head(all_mean_counts))

##definition of a pause from [Gillen et al 2021](https://doi.org/10.1186/s13059-021-02494-w)
#### a pause site in each condition was defined as a position with a RPF peak height ten times greater than the average RPF peak on the mRNA####

#calculate the average RPF occupancy / mRNA in cdss - control
all_mean_counts %>%
  filter(condition %in% c("WT", "KO")) %>%
  distinct() %>%
  mutate(codon = ceiling(cds_position / 3)) %>% #calculate codons
  arrange(transcript, condition, replicate, cds_position) %>% 
  group_by(transcript, condition, codon) %>%
  summarise(codon_normalised_cpm = mean(normalised_CPM), codon_seq = paste(Nucleotide, collapse="")) %>% #average the normalised cpm over the 3 nt of the codons and save the codon sequences
  ungroup() %>%
  distinct() %>%
  mutate(codon_seq = str_sub(codon_seq, 1, 3)) %>% #codon sequences are repeated 3 times (because we summarised 3 replicates just before) so just keep the 1st 3 letters (1st repeat)
  pivot_wider(names_from = condition, values_from = codon_normalised_cpm, names_glue = "{condition}_{.value}") -> selected_counts #create a column for each condition

print(head(selected_counts))
print(summary(selected_counts))

#calculate the average RPF peak size / mRNA in cdss - WT
selected_counts %>%
  filter(WT_codon_normalised_cpm > 0) %>%
  group_by(transcript) %>%
  summarise(WT_mean_height = mean(WT_codon_normalised_cpm)) -> wt_average_peak_height_codons

#calculate the average RPF peak size / mRNA in cdss - KO
selected_counts %>%
  filter(KO_codon_normalised_cpm > 0) %>%
  group_by(transcript) %>%
  summarise(KO_mean_height = mean(KO_codon_normalised_cpm)) -> ko_average_peak_height_codons

#merge with main file
selected_counts %>%
  left_join(wt_average_peak_height_codons, by = "transcript") %>%
  left_join(ko_average_peak_height_codons, by = "transcript") %>%
#identify pauses
  mutate(WT_pause = factor(case_when(!is.na(WT_mean_height) & WT_codon_normalised_cpm >= 10*WT_mean_height ~ 'pause',
                                         TRUE ~ 'no_pause')),
         KO_pause = factor(case_when(!is.na(KO_mean_height) & KO_codon_normalised_cpm >= 10*KO_mean_height ~ 'pause',
                                         TRUE ~ 'no_pause'))) -> all_mean_counts_pause_codon

#print the number of pauses found
print("NPM WT")
length(unique(all_mean_counts_pause_codon$transcript[all_mean_counts_pause_codon$WT_pause == "pause"]))
print("NPM KO")
length(unique(all_mean_counts_pause_codon$transcript[all_mean_counts_pause_codon$KO_pause == "pause"]))


all_mean_counts_pause_codon %>%
  inner_join(most_abundant_transcripts, by = "transcript") %>%
  filter(WT_pause == "pause" | KO_pause == "pause") -> all_mean_counts_pause_codon_filtered

####b. pause sites classification:####
#?sustained?: they are present but unaltered with treatment;
#"resolved? in that the reduction in RPF peak height is ten times greater than the average delta decrease across the mRNA 
#?induced? in that the increase in RPF peak height is ten times greater than the average delta increase

#calculate the average change in RPFs for each transcript
selected_counts %>%
#calculate delta(normalised_cpm)
  mutate(delta_cpm = KO_codon_normalised_cpm - WT_codon_normalised_cpm) %>%
#calculate the average change upon treatment across whole mRNAs in cdss 
  group_by(transcript) %>%
  mutate(average_delta = mean(delta_cpm)) %>%
  ungroup()  -> all_mean_counts_codons_delta

#calculate the average increase in peak heights for up genes 
all_mean_counts_codons_delta %>%
  filter(average_delta > 0) %>%
  filter(delta_cpm > 0) %>%
  group_by(transcript) %>%
  summarise(average_change = mean(delta_cpm)) -> increase

#calculate the average decreased in peak height for down genes 
all_mean_counts_codons_delta %>%
  filter(average_delta < 0) %>%
  filter(delta_cpm < 0) %>%
  group_by(transcript) %>%
  summarise(average_change = mean(delta_cpm)) -> decrease

#concatenate and merge with delta data
all_mean_counts_codons_delta %>%
  inner_join(most_abundant_transcripts, by = "transcript") %>%
  inner_join(rbind(increase, decrease), by = "transcript") %>%
  inner_join(all_mean_counts_pause_codon_filtered %>% select(all_of(c("transcript", "codon", "WT_pause", "KO_pause"))), by = c("transcript", "codon")) %>%
#define induced, resolved and sustained pauses
  mutate(pause_type = factor(case_when(average_delta > 0 & delta_cpm > 10*average_change ~ "induced", 
                                       average_delta < 0 & delta_cpm > -10*average_change ~ "induced", 
                                       average_delta < 0 & delta_cpm < 10*average_change ~ "resolved", 
                                       average_delta > 0 & delta_cpm < - 10*average_change ~ "resolved", 
                                       WT_pause == "pause" & KO_pause == "pause" ~ "sustained",
                                       TRUE ~ "ns")))  -> pause_type

#extract the sequence of the codons in the E and and A sites
pause_type %>%
  rename(p_site = codon, p_site_seq = codon_seq) %>%
  mutate(e_site = p_site - 1, 
         a_site = p_site + 1) %>%
  left_join(selected_counts %>% rename(e_site = codon, e_site_seq = codon_seq) %>% select(c("transcript", "e_site", "e_site_seq")), by = c("transcript", "e_site")) %>%
  left_join(selected_counts %>% rename(a_site = codon, a_site_seq = codon_seq) %>% select(c("transcript", "a_site", "a_site_seq")), by = c("transcript", "a_site")) %>%
  select(-c("e_site", "a_site")) -> pause_type

#export data
write_csv(pause_type, file.path(parent_dir, "NPM_pause_analysis.csv"))


