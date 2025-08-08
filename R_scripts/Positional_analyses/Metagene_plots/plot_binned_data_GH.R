# load libraries & set conflicts ------
library(tidyverse)
library(grid)
library(gridExtra)
library(viridis)
library(conflicted)

conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::select)
conflicts_prefer(purrr::reduce)

#read in common variables & shared files----
source("common_variables.R")
signal_P_data <- read_csv(file = paste0(machine_dir, "/R11/bioinformatics_resources/FASTAs/mouse/GENCODE/vM27/transcript_info/SignalP_6_prediction_results.csv"))
region_lengths <- read_csv(file = paste0(machine_dir, "/R11/bioinformatics_resources/FASTAs/mouse/GENCODE/vM27/transcript_info/gencode.vM27.pc_transcripts_region_lengths.csv"), 
                           col_names = c("transcript", "UTR5_len", "CDS_len", "UTR3_len"))
most_abundant_transcripts <- read_csv(file = file.path(parent_dir, "Analysis/most_abundant_transcripts/most_abundant_transcripts_IDs.csv"))

#set what you have called your control and treated samples. This can be a vector of strings if more than one treatment has been used.
control <- "WT"
treatment <- "KO"
TMT <- "NPM KO"
Tissue <- "APC-KRAS"

#read in functions----
source("binning_RiboSeq_functions.R")

#create themes----
my_theme <- theme_bw()+
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title = element_blank())

UTR5_theme <- my_theme+
  theme(legend.position="none",
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 18),
        axis.text.x = element_blank())

CDS_theme <- my_theme+
  theme(legend.position="none",
        axis.ticks = element_blank(),
        axis.text = element_blank())

UTR3_theme <- my_theme+
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.text = element_text(size = 18),
        legend.title = element_blank())

#read in data----
load(file = file.path(parent_dir, "Counts_files/R_objects/binned_list.Rdata"))
summary(binned_list[[1]])
head(binned_list[[1]])

load(file = file.path(parent_dir, "Counts_files/R_objects/single_nt_list.Rdata"))
summary(single_nt_list[[1]])
head(single_nt_list[[1]])

#all transcripts----
#summarise
#summarise within each sample
summarised_binned_list <- lapply(binned_list, summarise_data, value = "binned_cpm", grouping = "bin")
summary(summarised_binned_list[[1]])
print(summarised_binned_list[[1]])

#summarise within each condition (across replicates)
do.call("rbind", summarised_binned_list) %>%
  group_by(grouping, condition, region) %>%
  summarise(average_counts = mean(mean_counts),
            sd_counts = sd(mean_counts)) %>%
  ungroup() -> summarised_binned
summary(summarised_binned)
print(summarised_binned)

#plot lines
binned_line_plots <- plot_binned_lines(df = summarised_binned, SD = T, control = control, treatment = treatment, colors = c(WT_color, KO_color))

png(filename = file.path(parent_dir, "plots/binned_plots/all_transcripts/all transcripts binned lines.png"), width = 1000, height = 310)
grid.arrange(binned_line_plots[[1]], binned_line_plots[[2]], binned_line_plots[[3]], nrow = 1, widths = c(1,2,1.5),
             top = textGrob(paste(Tissue, TMT, "all transcripts binned"),gp=gpar(fontsize=16,font=2)))
dev.off()

binned_line_plots_all_replicates <- plot_binned_all_replicates(summarised_binned_list, control = control, treatment = treatment, colors = c(WT_color, KO_color))

png(filename = file.path(parent_dir, "plots/binned_plots/all_transcripts/all transcripts binned lines all replicates.png"), width = 1000, height = 310)
grid.arrange(binned_line_plots_all_replicates[[1]], binned_line_plots_all_replicates[[2]], binned_line_plots_all_replicates[[3]], nrow = 1, widths = c(1,2,1.5),
             top = textGrob(paste(Tissue, TMT, "all transcripts binned"),gp=gpar(fontsize=16,font=2)))
dev.off()

#calculate and plot delta
binned_delta_data <- calculate_binned_delta(binned_list, value = "binned_cpm", control = control, treatment = treatment, paired_data = F)
binned_delta_plots <- plot_binned_delta(binned_delta_data)

png(filename = file.path(parent_dir, "plots/binned_plots/all_transcripts/all transcripts binned delta.png"), width = 1000, height = 210)
grid.arrange(binned_delta_plots[[1]], binned_delta_plots[[2]], binned_delta_plots[[3]],
             nrow = 1, widths = c(1,2,1),
             top = textGrob(paste(Tissue, TMT, "all transcripts binned"),gp=gpar(fontsize=16,font=2)))
dev.off()

#positional----
#normalise within each transcript
positional_list <- lapply(binned_list, calculate_positional_counts)

#summarise within each sample
summarised_positional_list <- lapply(positional_list, summarise_data, value = "positional_counts", grouping = "bin")

#summarise within each condition (across replicates)
do.call("rbind", summarised_positional_list) %>%
  group_by(grouping, condition) %>%
  summarise(average_counts = mean(mean_counts),
            sd_counts = sd(mean_counts)) %>%
  ungroup() -> summarised_positional

#plot
positional_line_plots <- plot_positional_lines(df = summarised_positional, SD = T, control = control, treatment = treatment, colors = c(WT_color, KO_color))

png(filename = file.path(parent_dir, "plots/binned_plots/all_transcripts/all transcripts binned positional lines.png"), width = 500, height = 200)
print(positional_line_plots)
dev.off()

#calculate and plot delta
binned_positional_delta <- calculate_positional_delta(positional_list, control = control, treatment = treatment, paired_data = F)
positional_binned_delta_plots <- plot_positional_delta(binned_positional_delta)

png(filename = file.path(parent_dir, "plots/binned_plots/all_transcripts/all transcripts binned positional delta.png"), width = 500, height = 200)
print(positional_binned_delta_plots)
dev.off()

#single nt----
#summarise
summarised_single_nt_list <- lapply(single_nt_list, summarise_data, value = "single_nt_cpm", grouping = "window")
summary(summarised_single_nt_list[[1]])
print(summarised_single_nt_list[[1]])

do.call("rbind", summarised_single_nt_list) %>%
  group_by(grouping, condition, region) %>%
  summarise(average_counts = mean(mean_counts),
            sd_counts = sd(mean_counts)) %>%
  ungroup() -> summarised_single_nt
summary(summarised_single_nt)
print(summarised_single_nt)

#plot
single_nt_line_plots <- plot_single_nt_lines(summarised_single_nt, SD=T, plot_ends=F, control = control, treatment = treatment, colors = c(WT_color, KO_color))

png(filename = file.path(parent_dir, "plots/binned_plots/all_transcripts/all transcripts single nt lines.png"), width = 1300, height = 310)
grid.arrange(single_nt_line_plots[[1]], single_nt_line_plots[[2]], single_nt_line_plots[[3]], single_nt_line_plots[[4]], nrow = 1, widths = c(1,2,2,1.5),
             top = textGrob(paste(Tissue, TMT, "all transcripts single nt"),gp=gpar(fontsize=16,font=2)))
dev.off()

#calculate and plot delta
single_nt_delta_data <- calculate_single_nt_delta(single_nt_list, value = "single_nt_cpm", control = control, treatment = treatment, paired_data = F)
single_nt_delta_plots <- plot_single_nt_delta(single_nt_delta_data, SD = T)

png(filename = file.path(parent_dir, "plots/binned_plots/all_transcripts/all transcripts single nt delta.png"), width = 1300, height = 210)
grid.arrange(single_nt_delta_plots[[1]], single_nt_delta_plots[[2]], single_nt_delta_plots[[3]], single_nt_delta_plots[[4]], nrow = 1, widths = c(1,2,2,1),
             top = textGrob(paste(Tissue, TMT, "all transcripts single nt"),gp=gpar(fontsize=16,font=2)))
dev.off()

# By Polarity score -----
df <- read_tsv(file = file.path(parent_dir,"RiboMiner/NPM_KO_polarity_dataframe.txt"))

df %>%
  rename(transcript = 1) %>%
  inner_join(most_abundant_transcripts, by = "transcript") %>%
  relocate(c("gene", "gene_sym"), .after = "transcript") %>%
  rowwise() %>%
  mutate(ave_WT = base::mean(c_across(c(4:7)), na.rm = TRUE),
         ave_KO = base::mean(c_across(c(8:11)), na.rm = TRUE),
         Delta_Polarity = ave_KO - ave_WT) -> df_delta

df_delta <- as.data.frame(df_delta)
df_delta <- as_tibble(df_delta)

## by Delta polarity score ----

df_delta %>%
  arrange(Delta_Polarity) %>%
  slice_head(prop = 0.05) %>%
  pull(transcript) -> PS_bottom_5

df_delta %>%
  arrange(Delta_Polarity) %>%
  slice_tail(prop = 0.05) %>%
  pull(transcript) -> PS_top_5

df_delta %>%
  arrange(Delta_Polarity) %>%
  slice_head(prop = 0.01) %>%
  pull(transcript) -> PS_bottom_1

df_delta %>%
  arrange(Delta_Polarity) %>%
  slice_tail(prop = 0.01) %>%
  pull(transcript) -> PS_top_1

df_delta %>%
  arrange(Delta_Polarity) %>%
  slice_head(prop = 0.1) %>%
  pull(transcript) -> PS_bottom_10

df_delta %>%
  arrange(Delta_Polarity) %>%
  slice_tail(prop = 0.1) %>%
  pull(transcript) -> PS_top_10

#plot binned
plot_subset(IDs = PS_bottom_5, subset = "Delta_bottom5perc", sub_dir = "Polarity_score",
            control = control, treatment = treatment, colors = c(WT_color, KO_color),
            binned_value = "binned_normalised_cpm", single_nt_value = "single_nt_normalised_cpm",
            plot_binned = T, plot_single_nt = F, plot_positional = T,
            plot_replicates = T, plot_delta = T, SD = T, paired_data = F)

plot_subset(IDs = PS_top_5, subset = "Delta_top5perc", sub_dir = "Polarity_score",
            control = control, treatment = treatment, colors = c(WT_color, KO_color),
            binned_value = "binned_normalised_cpm", single_nt_value = "single_nt_normalised_cpm",
            plot_binned = T, plot_single_nt = F, plot_positional = T,
            plot_replicates = T, plot_delta = T, SD = T, paired_data = F)

plot_subset(IDs = PS_bottom_1, subset = "Delta_bottom1perc", sub_dir = "Polarity_score",
            control = control, treatment = treatment, colors = c(WT_color, KO_color),
            binned_value = "binned_normalised_cpm", single_nt_value = "single_nt_normalised_cpm",
            plot_binned = T, plot_single_nt = F, plot_positional = T,
            plot_replicates = T, plot_delta = T, SD = T, paired_data = F)

plot_subset(IDs = PS_top_1, subset = "Delta_top1perc", sub_dir = "Polarity_score",
            control = control, treatment = treatment, colors = c(WT_color, KO_color),
            binned_value = "binned_normalised_cpm", single_nt_value = "single_nt_normalised_cpm",
            plot_binned = T, plot_single_nt = F, plot_positional = T,
            plot_replicates = T, plot_delta = T, SD = T, paired_data = F)

plot_subset(IDs = PS_bottom_10, subset = "Delta_bottom10perc", sub_dir = "Polarity_score",
            control = control, treatment = treatment, colors = c(WT_color, KO_color),
            binned_value = "binned_normalised_cpm", single_nt_value = "single_nt_normalised_cpm",
            plot_binned = T, plot_single_nt = F, plot_positional = T,
            plot_replicates = T, plot_delta = T, SD = T, paired_data = F)

plot_subset(IDs = PS_top_10, subset = "Delta_top10perc", sub_dir = "Polarity_score",
            control = control, treatment = treatment, colors = c(WT_color, KO_color),
            binned_value = "binned_normalised_cpm", single_nt_value = "single_nt_normalised_cpm",
            plot_binned = T, plot_single_nt = F, plot_positional = T,
            plot_replicates = T, plot_delta = T, SD = T, paired_data = F)

## By Polarity score per condition

df_delta %>%
  arrange(ave_WT) %>%
  slice_head(prop = 0.05) %>%
  pull(transcript) -> bottom_5_WT

df_delta %>%
  arrange(ave_KO) %>%
  slice_head(prop = 0.05) %>%
  pull(transcript) -> bottom_5_KO

df_delta %>%
  arrange(ave_WT) %>%
  slice_tail(prop = 0.05) %>%
  pull(transcript) -> top_5_WT

df_delta %>%
  arrange(ave_KO) %>%
  slice_tail(prop = 0.05) %>%
  pull(transcript) -> top_5_KO

plot_subset(IDs = bottom_5_WT, subset = "Bottom5_WT", sub_dir = "Polarity_score",
            control = control, treatment = treatment, colors = c(WT_color, KO_color),
            binned_value = "binned_normalised_cpm", single_nt_value = "single_nt_normalised_cpm",
            plot_binned = T, plot_single_nt = F, plot_positional = T,
            plot_replicates = T, plot_delta = T, SD = T, paired_data = F)

plot_subset(IDs = bottom_5_KO, subset = "Bottom5_KO", sub_dir = "Polarity_score",
            control = control, treatment = treatment, colors = c(WT_color, KO_color),
            binned_value = "binned_normalised_cpm", single_nt_value = "single_nt_normalised_cpm",
            plot_binned = T, plot_single_nt = F, plot_positional = T,
            plot_replicates = T, plot_delta = T, SD = T, paired_data = F)

plot_subset(IDs = top_5_WT, subset = "Top5_WT", sub_dir = "Polarity_score",
            control = control, treatment = treatment, colors = c(WT_color, KO_color),
            binned_value = "binned_normalised_cpm", single_nt_value = "single_nt_normalised_cpm",
            plot_binned = T, plot_single_nt = F, plot_positional = T,
            plot_replicates = T, plot_delta = T, SD = T, paired_data = F)

plot_subset(IDs = top_5_KO, subset = "Top5_KO", sub_dir = "Polarity_score",
            control = control, treatment = treatment, colors = c(WT_color, KO_color),
            binned_value = "binned_normalised_cpm", single_nt_value = "single_nt_normalised_cpm",
            plot_binned = T, plot_single_nt = F, plot_positional = T,
            plot_replicates = T, plot_delta = T, SD = T, paired_data = F)

# GSEA pathways----
library(fgsea)

## read in pathways - perso folder ----
source("read_mouse_GSEA_pathways.R")

# function going through again the groups of genes in the leading edges of proteomics and RPFs
plot_subsets_shared <- function(pathway) {
  
  # get info out from the leading edge merged data
  LE %>% filter(pathway == !!pathway) -> POI_row
  
  prot_genes <- if(nchar(POI_row$leadingEdge_Prot) > 0) {
    str_split(POI_row$leadingEdge_Prot, "\\s+")[[1]]
  } else {
    character(0)
  }
  
  rpf_genes <- if(nchar(POI_row$leadingEdge_RPFs) > 0) {
    str_split(POI_row$leadingEdge_RPFs, "\\s+")[[1]]
  } else {
    character(0)
  }
  
  common_genes <- intersect(prot_genes, rpf_genes)
  only_prot <- setdiff(prot_genes, rpf_genes)
  only_rpf <- setdiff(rpf_genes, prot_genes)
  
  # Assuming gene_transcript_table has columns: gene and transcript
  common_transcripts <- most_abundant_transcripts$transcript[match(common_genes, most_abundant_transcripts$gene_sym)]
  only_prot_transcripts <- most_abundant_transcripts$transcript[match(only_prot, most_abundant_transcripts$gene_sym)]
  only_rpf_transcripts  <- most_abundant_transcripts$transcript[match(only_rpf,  most_abundant_transcripts$gene_sym)]
  
  plot_subset(IDs = common_transcripts, subset = paste0("common_", pathway), sub_dir = "LeadingEdges",
              control = control, treatment = treatment, colors = c(WT_color, KO_color),
              binned_value = "binned_cpm", single_nt_value = "single_nt_normalised_cpm",
              plot_binned = T, plot_single_nt = F, plot_positional = T,
              plot_replicates = T, plot_delta = T, SD = T, paired_data = F)
  
  plot_subset(IDs = only_prot_transcripts, subset = paste0("OnlyProt_", pathway), sub_dir = "LeadingEdges",
              control = control, treatment = treatment, colors = c(WT_color, KO_color),
              binned_value = "binned_cpm", single_nt_value = "single_nt_normalised_cpm",
              plot_binned = T, plot_single_nt = F, plot_positional = T,
              plot_replicates = T, plot_delta = T, SD = T, paired_data = F)
  
  plot_subset(IDs = only_rpf_transcripts, subset = paste0("OnlyRPFs_", pathway), sub_dir = "LeadingEdges",
              control = control, treatment = treatment, colors = c(WT_color, KO_color),
              binned_value = "binned_cpm", single_nt_value = "single_nt_normalised_cpm",
              plot_binned = T, plot_single_nt = F, plot_positional = T,
              plot_replicates = T, plot_delta = T, SD = T, paired_data = F)
  
}

# read in leading edges HALLMARKS ---------
LE_RPFs <- read_tsv(file = file.path(parent_dir, "plots/fgsea/scatters/AK_Npm1KO_Hallmark_RPF_BC_LeadEdge.tsv")) 
LE_Protein <- read_tsv(file = file.path(parent_dir, "plots/fgsea/Proteomics/scatters/AK_Npm1KO_Hallmark_LeadEdge.tsv"))

HM <- data.frame(pathway = names(pathways.hallmark),
                 count = sapply(pathways.hallmark, length),
                 stringsAsFactors = FALSE)
HM %>% remove_rownames() -> HM

# as of 20250325 - this analysis done on all shared pathways and not only on the significant ones
LE_Protein %>% 
  left_join(LE_RPFs, by = join_by(pathway), suffix = c("_Prot", "_RPFs")) %>% 
  left_join(HM)-> LE

LE %>% pull(pathway) -> POI

lapply(POI, function(pathway) {
  tryCatch({
    plot_subsets_shared(pathway)
  }, error = function(e) {
    warning(paste("Plotting for pathway", pathway, "failed with error:", e$message))
  })
})

# read in leading edges GOBP ---------
LE_RPFs <- read_tsv(file = file.path(parent_dir, "plots/fgsea/scatters/AK_Npm1KO_GOBP_RPF_BC_LeadEdge.tsv")) 
LE_Protein <- read_tsv(file = file.path(parent_dir, "plots/fgsea/Proteomics/scatters/AK_Npm1KO_BP_LeadEdge.tsv"))

GOBP <- data.frame(pathway = names(pathways.bio_processes),
                   count = sapply(pathways.bio_processes, length),
                   stringsAsFactors = FALSE)
GOBP %>% remove_rownames() -> GOBP

# as of 20250325 - this analysis done on all shared pathways and not only on the significant ones
LE_Protein %>% 
  left_join(LE_RPFs, by = join_by(pathway), suffix = c("_Prot", "_RPFs")) %>% 
  left_join(GOBP) -> LE

# for GOBP there are so many pathways that it makes no sense to plot all
# plotting is done only for the two pathways of interest, which are ER stress and ISR

LE %>%
  filter(pathway == "GOBP_RESPONSE_TO_ENDOPLASMIC_RETICULUM_STRESS" | 
           pathway == "GOBP_INTEGRATED_STRESS_RESPONSE_SIGNALING") %>% 
  pull(pathway) -> POI

lapply(POI, function(pathway) {
  tryCatch({
    plot_subsets_shared(pathway)
  }, error = function(e) {
    warning(paste("Plotting for pathway", pathway, "failed with error:", e$message))
  })
})
