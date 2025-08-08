# We saw in the GSEAs that some of the Hallmark pathways had an opposite direction of enrichments (see Fig4)
# We therefore wanted to check if indeed the genes belonging to the leading edges (ie the genes driving the enrichment)
# are more or less abundant in one or the other omics

# The bottom of the script adds also the total cytoplamic RNA data scaled
# Was primarily used to check for the KD and for the p53 data in the paper

# needs:
# Required data for plots
# RPFs DE analysis
# RPFs normalised reads
# Leading edges for Hallmarks of Proteomics and RPFs
# Proteomics - needs more parsing than other tables

# In this script the definition of Shared / Only Protein / Only RPFs referes to whether the gene is part of the Leading Edge for the specific pathway

# load libraries & set conflicts ------
library(tidyverse)
library(ggrepel)
library(readxl)
library(conflicted)
library(ggh4x)

conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::select)
conflicts_prefer(purrr::reduce)

# make sure you are in the R_script folder then call
source("common_variables.R")

plot_dir <- paste0(parent_dir, "/plots/fgsea/Proteomics_RPF_comparison")

control <- "WT"
treatment <- "KO"

Tissue <- "APC-KRAS"

# themes & functions ----
mytheme <- theme_classic()+
  theme(plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14))

# Plotting functions -----
plot_vioilns_groups <- function(df, pathway) {
  
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
  
  POI_row %>% pull(NES_Prot) -> N_P
  POI_row %>% pull(NES_RPFs) -> N_R
  POI_row %>% pull(padj_Prot) -> padj_P
  POI_row %>% pull(padj_RPFs) -> padj_R
  
  POI_row %>% pull(count) -> GIP # genes in pathway
  POI_row %>% pull(size_RPFs) -> DR # detected ion riboseq 
  POI_row %>% pull(size_Prot) -> DP # detected in proteomics
  
  df %>%
    mutate(gene_group = case_when(gene_sym %in% common_genes ~ "Shared",
                                  gene_sym %in% only_prot ~ "Only Protein",
                                  gene_sym %in% only_rpf ~ "Only RPFs",
                                  TRUE ~ "none"),
           alpha_score = ifelse(gene_group %in% c("Shared", "Only Protein", "Only RPFs"), 1, 0.01)) %>%
    mutate(gene_group = factor(gene_group, levels = c("Shared", "Only Protein", "Only RPFs", "none"))) %>%
    arrange(desc(gene_group)) -> df

  df %>% 
    filter(gene_group != "none") %>% 
    distinct(gene_sym, gene_group) %>% 
    count(gene_group) -> group_counts
  
  setNames(paste0(group_counts$gene_group, " (n=", group_counts$n, ")"),
    group_counts$gene_group) -> labeller_vec
  
  df %>% 
    filter(gene_group != "none") %>% 
    pivot_longer(cols = -c(gene_sym, gene_group, alpha_score),
                 names_to = c("assay", "condition", "replicate"),
                 names_pattern = "(Protein|RPFs)_([A-Z]{2})_(\\d+)") %>% 
    mutate(group = case_when(
      assay == "Protein" & condition == "WT" ~ "WT Protein",
      assay == "Protein" & condition == "KO" ~ "KO Protein",
      assay == "RPFs"  & condition == "WT" ~ "WT RPFs",
      assay == "RPFs"  & condition == "KO" ~ "KO RPFs")) %>% 
    mutate(group = factor(group, levels = c("WT RPFs", "KO RPFs", "WT Protein", "KO Protein"))) %>%
    filter(!is.na(group)) %>% 
    ggplot(aes(x = group, y = value, fill = group)) +
    geom_violin(trim = FALSE, linewidth = 0.6) +
    geom_boxplot(position = position_dodge(0.8), width = 0.2, outlier.shape = NA, linewidth = 0.4) +
    facet_wrap(~ gene_group, scales = "free_y", labeller = as_labeller(labeller_vec)) +
    scale_fill_manual(values = c("grey", "firebrick3", "darkgrey", "darkred")) +
    labs(
      x = "Measurement Group",
      y = "Normalized Quantification",
      title = paste(paste(Tissue, "Npm1 KO"), str_replace_all(pathway, "_", " "), paste("Genes in pathway:", GIP), sep = "\n"),
      subtitle = paste("NES Protein", round(N_P, 2),
                       "\nadj pvalue Protein", signif(padj_P, 2),
                       "\nDetected Proteins:", DP,
                       "\nNES RPFs", round(N_R, 2),
                       "\nadj pvalue RPFs", signif(padj_R, 2),
                       "\nDetected in Riboseq:", DR)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_text(size = 12),
          axis.title.y = element_text(size = 12),
          strip.text = element_text(size = 14),
          legend.position = "none",
    ) -> violin
  
  pdf(file = paste0(plot_dir, "/violins/", paste(pathway, "RPF_Prot_violin_plot_AKonly.pdf", sep = "_")), width = 10, height = 7)
  print(violin)
  dev.off()
  
  return(violin)
}

plot_vioilns_shared <- function(df, pathway) {
  
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
  
  POI_row %>% pull(NES_Prot) -> N_P
  POI_row %>% pull(NES_RPFs) -> N_R
  POI_row %>% pull(padj_Prot) -> padj_P
  POI_row %>% pull(padj_RPFs) -> padj_R
  
  POI_row %>% pull(count) -> GIP # genes in pathway
  POI_row %>% pull(size_RPFs) -> DR # detected ion riboseq
  POI_row %>% pull(size_Prot) -> DP # detected in proteomics
  
  df %>%
    mutate(gene_group = case_when(gene_sym %in% common_genes ~ "Shared",
                                  TRUE ~ "none"),
           alpha_score = ifelse(gene_group %in% "Shared", 1, 0.01)) %>%
    mutate(gene_group = factor(gene_group)) %>%
    arrange(desc(gene_group)) -> df
  
  df %>%
    filter(gene_group == "Shared") %>% 
    tally() %>% 
    pull(n) -> GS
  
  df %>% 
    filter(gene_group != "none") %>% 
    pivot_longer(cols = -c(gene_sym, gene_group, alpha_score),
                 names_to = c("assay", "condition", "replicate"),
                 names_pattern = "(Protein|RPFs)_([A-Z]{2})_(\\d+)") %>% 
    mutate(group = case_when(
      assay == "Protein" & condition == "WT" ~ "WT Protein",
      assay == "Protein" & condition == "KO" ~ "KO Protein",
      assay == "RPFs"  & condition == "WT" ~ "WT RPFs",
      assay == "RPFs"  & condition == "KO" ~ "KO RPFs")) %>% 
    mutate(group = factor(group, levels = c("WT RPFs", "KO RPFs", "WT Protein", "KO Protein"))) %>%
    filter(!is.na(group)) %>% 
    ggplot(aes(x = group, y = value, fill = group)) +
    geom_violin(trim = FALSE, linewidth = 0.6) +
    geom_boxplot(position = position_dodge(0.8), width = 0.2, outlier.shape = NA, linewidth = 0.4) +
    scale_fill_manual(values = c("grey", "firebrick3", "darkgrey", "darkred")) +
    labs(y = "Normalized Quantification",
         title = paste(paste0(Tissue, " Npm1 KO\n", str_wrap(str_replace_all(pathway, "_", " "), width = 30)),
                          paste("Genes in pathway:", GIP), 
                          paste("Shared genes in Leading Edges:", GS), sep = "\n"),
         subtitle = paste("NES Protein", round(N_P, 2),
                       "\nadj pvalue Protein", signif(padj_P, 2),
                       "\nDetected Proteins:", DP,
                       "\nNES RPFs", round(N_R, 2),
                       "\nadj pvalue RPFs", signif(padj_R, 2),
                       "\nDetected in Riboseq:", DR)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_text(size = 12),
          axis.title.y = element_text(size = 12),
          legend.position = "none",
    ) -> violin
  
  pdf(file = paste0(plot_dir, "/violins/", paste(pathway, "RPF_Prot_violin_plot_onlyShared_AKonly.pdf", sep = "_")), width = 5.5, height = 6.5)
  print(violin)
  dev.off()
  
  return(violin)
}

# read in Riboseq - merged results -------
DE <- read_csv(file = file.path(parent_dir, "Analysis/DESeq2_output", "merged_DESeq2_AK_NPM_KO_batchCorrected_unshrunk.csv")) # this contains also totals and TE

# read in RPFs normalised counts -----
RPFs_nc <- read_csv(file = file.path(parent_dir, "Analysis/DESeq2_output", "RPFs_APC-KRAS_KO_normalised_counts.csv"))

# read in totals normalised counts -----
Tots_nc <- read_csv(file = file.path(parent_dir, "Analysis/DESeq2_output", "totals_APC-KRAS_KO_normalised_counts.csv"))

# read in proteomics results from excel file ------
TMT <- read_excel(paste0(machine_dir, "/CGIACOME/AAJ_NPM/20230816_TMT_proteins/20231030_results/CHG_R11_160823_updated.xlsx"), sheet = "Proteome_Main")

TMT %>% 
  rename("Gene_ID" = 1,
         "Protein_ID" = 2,
         "Protein_names" = 3,
         "gene_sym" = 4,
         "Ttest_AK_WToverKO" = 5,
         "Ttest_APC_WToverKO" = 6,
         "Ttest_sign" = 7,
         "pval_AK_WToverKO" = 8,
         "pval_APC_WToverKO" = 9,
         "Sign_AK" = 10,
         "Sign_APC" = 11,
         "Anova_sig" = 12,
         "Anova_Cluster" = 13,
         "AK_WT_1" =  14,
         "AK_WT_2" = 15,
         "AK_WT_3" = 16,
         "AK_WT_4" = 17,
         "AK_KO_1" = 18,
         "AK_KO_2" = 19,
         "AK_KO_3" = 20,
         "AK_KO_4" = 21,
         "APC_WT_1" = 22,
         "APC_WT_2" = 23,
         "APC_WT_3" = 24,
         "APC_WT_4" = 25,
         "APC_KO_1" = 26,
         "APC_KO_2" = 27,
         "APC_KO_3" = 28,
         "APC_KO_4" = 29) -> TMT

TMT %>% 
  select(c(2, 4, 14:21)) %>% 
  filter(!is.na(gene_sym)) -> TMT_for_HM

TMT_for_HM %>% 
  pivot_longer(cols = c(3:10)) %>% 
  ggplot(aes(x = value, color = name)) +
  geom_density()

TMT_for_HM %>% 
  group_by(Protein_ID) %>% 
  filter(n() >1) %>% 
  summarise(n = n()) -> duplicate_prot # check that there is no duplicate proteins by ID

TMT_for_HM %>% 
  select(c(-2)) %>% 
  column_to_rownames("Protein_ID") -> to_transform

Scaled_for_HM <- t(scale(t(to_transform)))

Scaled_for_HM %>% 
  as_tibble(rownames = NA) %>% 
  rownames_to_column("Protein_ID") %>% 
  inner_join(TMT_for_HM[c(1,2)]) %>% 
  relocate(gene_sym, .after = Protein_ID) -> Scaled_data

summary(Scaled_data)

Scaled_data %>% 
  filter(gene_sym == "Npm1") # cross check KO

Scaled_data %>% 
  rowwise() %>% 
  mutate(ave_AK_WT = mean(c(AK_WT_1, AK_WT_2, AK_WT_3, AK_WT_4)),
         ave_AK_KO = mean(c(AK_KO_1, AK_KO_2, AK_KO_3, AK_KO_4))) %>% 
  relocate(c(11,12), .after = gene_sym) %>% 
  as_tibble() -> Scaled_data_final

# need to pivot longer sep where the protein group still show more than one protein to do a full match
# caveat: some proteins will have different effects at RPFs but the exact same effect on proteomics if they were in the same protein group
Scaled_data_final %>%
  left_join(TMT[c(2,5:21)], by = join_by(Protein_ID), suffix = c("_scaled", "_orig")) %>%
  separate_longer_delim(gene_sym, delim = ";") -> merged_longer

merged_longer %>% 
  rowwise() %>% 
  mutate(ave_AK_WT_orig = mean(c(AK_WT_1_orig, AK_WT_2_orig, AK_WT_3_orig, AK_WT_4_orig)),
         ave_AK_KO_orig = mean(c(AK_KO_1_orig, AK_KO_2_orig, AK_KO_3_orig, AK_KO_4_orig))) %>% 
  relocate(ave_AK_WT_orig, ave_AK_KO_orig, .after = gene_sym) %>% 
  group_by(gene_sym) %>% # get rid of duplicates and then check in quant_diff that none of them had differences than the removed
  mutate(
    avg_quant = mean(Ttest_AK_WToverKO, na.rm = TRUE),
    quant_diff = max(Ttest_AK_WToverKO, na.rm = TRUE) - min(Ttest_AK_WToverKO, na.rm = TRUE)
  ) %>%
  slice_max(order_by = ave_AK_WT_orig, n = 1, with_ties = FALSE) %>%
  ungroup() %>% 
  left_join(DE, multiple = "all", relationship = "many-to-many") -> merged_data

merged_data %>% 
  group_by(gene_sym) %>% 
  filter(n() >1) %>% 
  arrange(gene_sym) -> duplicate_gene_sym # A couple of gene symbols had different ENSEMBL identifiers

# let's cross check the KO
merged_data %>% 
  select(gene_sym, Ttest_AK_WToverKO, Ttest_APC_WToverKO, RPFs_log2FC, totals_log2FC) %>% 
  filter(gene_sym == "Npm1") # no panick, the reason ttest results are positive is just on how Proteomics CF did the ttest (with WT over KO, instead of KO over WT)

# Merge normalised counts after scaling them for the violin plots
RPFs_nc %>% 
  select(1:9) %>% 
  column_to_rownames("transcript") -> to_transform

Scaled_RPFs_nc <- t(scale(t(to_transform)))

Scaled_RPFs_nc %>% 
  as_tibble(rownames = NA) %>% 
  rownames_to_column("transcript") %>% 
  rowwise() %>% 
  mutate(ave_RPF_WT = mean(c(NPM_WT_1_RPFs, NPM_WT_2_RPFs, NPM_WT_3_RPFs, NPM_WT_4_RPFs)),
         ave_RPF_KO = mean(c(NPM_KO_1_RPFs, NPM_KO_2_RPFs, NPM_KO_3_RPFs, NPM_KO_4_RPFs))) %>% 
  inner_join(RPFs_nc[c("transcript", "gene", "gene_sym")], by = join_by("transcript")) %>% 
  relocate(gene, gene_sym, ave_RPF_WT, ave_RPF_KO, .after = transcript) %>% 
  as_tibble() %>% 
  arrange(desc(ave_RPF_WT)) -> Scaled_RPFs_nc

new_names <- c(
  paste0("Protein_WT_", 1:4),
  paste0("Protein_KO_", 1:4),
  paste0("RPFs_WT_", 1:4),
  paste0("RPFs_KO_", 1:4)
)

merged_data %>% 
  left_join(Scaled_RPFs_nc) %>% 
  select(c(2, 7:14, 46:53)) %>% 
  rename_with(~ new_names, .cols = 2:17)-> ready_4_violins

# let's cross check the KO
ready_4_violins %>% 
  filter(gene_sym == "Npm1")

# read in leading edges HALLMARKS ---------
LE_RPFs <- read_tsv(file = file.path(parent_dir, "plots/fgsea/scatters/AK_Npm1KO_Hallmark_RPF_BC_LeadEdge.tsv")) 
LE_Protein <- read_tsv(file = file.path(parent_dir, "plots/fgsea/Proteomics/scatters/AK_Npm1KO_Hallmark_LeadEdge.tsv"))

# fetch number of genes in each pathway - regardless of whether they are detected or not
source("read_mouse_GSEA_pathways.R")

HM <- data.frame(pathway = names(pathways.hallmark),
                 count = sapply(pathways.hallmark, length),
                 stringsAsFactors = FALSE)
HM %>% remove_rownames() -> HM

# as of 20250325 - this analysis done on all shared pathways and not only on the significant ones
LE_Protein %>% 
  left_join(LE_RPFs, by = join_by(pathway), suffix = c("_Prot", "_RPFs")) %>% 
  left_join(HM)-> LE

LE %>% pull(pathway) -> POI

# Now plot -----
# violins
lapply(POI, plot_vioilns_groups, df = ready_4_violins)

# violins only shared
lapply(POI, plot_vioilns_shared, df = ready_4_violins)

# read in leading edges GOBP ---------
LE_RPFs <- read_tsv(file = file.path(parent_dir, "plots/fgsea/scatters/AK_Npm1KO_GOBP_RPF_BC_LeadEdge.tsv")) 
LE_Protein <- read_tsv(file = file.path(parent_dir, "plots/fgsea/Proteomics/scatters/AK_Npm1KO_BP_LeadEdge.tsv"))

GOBP <- data.frame(pathway = names(pathways.bio_processes),
                 count = sapply(pathways.bio_processes, length),
                 stringsAsFactors = FALSE)
GOBP %>% remove_rownames() -> GOBP

LE_Protein %>% 
  left_join(LE_RPFs, by = join_by(pathway), suffix = c("_Prot", "_RPFs")) %>% 
  left_join(GOBP) -> LE

# for GOBP there are so many pathways that it makes no sense to plot all
# plotting is done only for the two pathways of interest, which are ER stress and ISR
LE %>%
  filter(pathway == "GOBP_RESPONSE_TO_ENDOPLASMIC_RETICULUM_STRESS" | 
           pathway == "GOBP_INTEGRATED_STRESS_RESPONSE_SIGNALING") %>% 
  pull(pathway) -> POI

# Now plot -----
# violins
lapply(POI, plot_vioilns_groups, df = ready_4_violins)

# violins only shared
lapply(POI, plot_vioilns_shared, df = ready_4_violins)

# merged totals normalised counts (after scaling) also to the ready 4 violins ------
# plot dot plots for individual genes
Tots_nc %>% 
  select(2:10) %>% 
  column_to_rownames("transcript") -> to_transform

Scaled_Tots_nc <- t(scale(t(to_transform)))

Scaled_Tots_nc %>% 
  as_tibble(rownames = NA) %>% 
  rownames_to_column("transcript") %>% 
  rowwise() %>%
  mutate(ave_Tot_WT = mean(c(NPM_WT_1_Totals, NPM_WT_2_Totals, NPM_WT_3_Totals, NPM_WT_4_Totals)),
         ave_Tot_KO = mean(c(NPM_KO_1_Totals, NPM_KO_2_Totals, NPM_KO_3_Totals, NPM_KO_4_Totals))) %>%
  inner_join(Tots_nc[c("transcript", "gene", "gene_sym")], by = join_by("transcript")) %>% 
  relocate(gene, gene_sym, ave_Tot_WT, ave_Tot_KO, .after = transcript) %>% 
  as_tibble() %>% 
  arrange(desc(ave_Tot_WT)) -> Scaled_Tots_nc

new_names <- c("Protein_WT_average", "Protein_KO_average",
               paste0("Protein_WT_", 1:4),
               paste0("Protein_KO_", 1:4),
               "RPFs_WT_average", "RPFs_KO_average",
               paste0("RPFs_WT_", 1:4),
               paste0("RPFs_KO_", 1:4),
               "Totals_WT_average", "Totals_KO_average",
               paste0("Totals_WT_", 1:4),
               paste0("Totals_KO_", 1:4)
               )

merged_data %>% 
  left_join(Scaled_RPFs_nc) %>% 
  left_join(Scaled_Tots_nc, by = join_by(transcript, gene)) %>% 
  select(c(2, 5:14, 43:53, 55:64)) %>% 
  relocate(transcript, .after = gene_sym.x) %>% 
  rename_with(~ new_names, .cols = 3:32) %>% 
  rename(gene_sym = gene_sym.x)-> ready_4_plots
  
# let's cross check the KO
ready_4_plots %>% 
  filter(gene_sym == "Npm1") 

# plot individual genes - for example check Npm1 KO across all omics
ready_4_plots %>% 
  select(-contains("_average")) %>% 
  filter(gene_sym == "Npm1") %>% # Choose here which gene you want to plot
  pivot_longer(
    cols = -c(gene_sym, transcript),
    names_to = c("assay", "condition", "replicate"),
    names_pattern = "(Protein|RPFs|Totals)_([A-Z]{2})_(\\d+)") %>% 
  mutate(assay = recode(assay, "Totals" = "Total RNA"),
         nested_label = paste(condition, assay, sep = "\n")) %>% 
  mutate(nested_label = factor(nested_label, 
                               levels = c("WT\nTotal RNA", "KO\nTotal RNA",
                                          "WT\nRPFs", "KO\nRPFs",
                                          "WT\nProtein", "KO\nProtein"))) %>%
  ggplot(aes(x = nested_label, y = value, fill = nested_label)) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.4) +
  geom_point(position = position_jitter(width = 0.15), size = 3) +
  scale_fill_manual(values = c("lightgrey", "red", "grey", "firebrick3", "darkgrey", "darkred")) +
  scale_x_discrete(labels = list( flat = rep(c("WT", "KO"), 3),
                                  nested = c("Total RNA", "RPFs", "Protein"),
                                  guide = guide_axis_nested())) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = "none",
  ) +
  labs(y = "Scaled Quantification", x = " sample group") -> plote

pdf(file = paste0(parent_dir, "/plots/", "Npm1_all_omics.pdf"), width = 4, height = 5) # Don't forget to change the title if you change gene name
print(plote)
dev.off()