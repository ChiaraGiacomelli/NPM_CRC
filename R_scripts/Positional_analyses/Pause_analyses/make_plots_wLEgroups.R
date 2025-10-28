# This script was originally developed by Dr.Pauline Herviou, Postdoctoral Scientist in the lab of Prof. Martin Bushell
# At the CRUK Scotland Institute, Glasgow

# The plot function from her script was then adapted to suit the leading edge data from the Gene Set Enrichment Analyses
# It looks at the presence and position of pause sites (sustained, resolved, or induced upon Npm1 KO) in the gene groups
# Where the gene groups are those in the Leading Edge of a pathway for both RPFs and Proteomics, or only for one of the two omics approaches

# load libraries & set conflicts ------
library(tidyverse)
library(scales)
library(gridExtra)
library(grid)
library(fgsea)
library(conflicted)

conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::select)
conflicts_prefer(purrr::reduce)

control <- "WT"
treatment <- "KO"

Tissue <- "APC-KRAS"

#read in common variables
source(paste0(parent_dir,"/R_scripts/common_variables.R"))

#read in functions----
source(paste0(parent_dir,"/R_scripts/binning_RiboSeq_functions.R"))
plot_dir <- paste0(parent_dir, "/plots/pauses/")

# import pause tables
# This is also shared here in the GitHub folder
pauses <- read_csv("NPM_pause_analysis_20250331.csv")

#import most abundant transcripts info
most_abundant_transcripts <- read_csv(file = "most_abundant_transcripts/most_abundant_transcripts_IDs.csv"))

#import region lengths table
region_lengths <- read_csv(file = "/FASTAs/mouse/GENCODE/vM27/transcript_info/gencode.vM27.pc_transcripts_region_lengths.csv",
                                         col_names = c("transcript", "UTR5_len", "CDS_len", "UTR3_len"))

# Codon properties tables are available here in the GitHub folder
codon_AA <- read_csv(file = "codon_amino_acid_pairs_and_properties.csv")
codon_anticod <- read_csv(file = "codon_anticodon_pairs.csv")

# Format table to have all info -----
pauses %>%
  mutate(across(c(e_site_seq, p_site_seq, a_site_seq), ~as.character(gsub("T", "U", .)))) %>% # make into "RNA" codons with U instead of T
  left_join(codon_AA, by = join_by(p_site_seq == codon)) %>% 
  left_join(codon_AA, by = join_by(a_site_seq == codon), suffix = c("_P", "_A")) %>% 
  left_join(codon_AA, by = join_by(e_site_seq == codon)) %>% 
  rename(AA_E = AA,
         Symbol_E = Symbol,
         Properties_E = Properties) %>% 
  mutate(
    E3 = str_sub(e_site_seq, -1),  # Extract last character of e_site_seq
    P3 = str_sub(p_site_seq, -1),  # Extract last character of p_site_seq
    A3 = str_sub(a_site_seq, -1)   # Extract last character of a_site_seq
  ) %>% 
  filter(!is.na(e_site_seq) & !is.na(a_site_seq) & !is.na(AA_A)) %>% # removes pause sites at the start codon and stop codon
  filter(pause_type != "ns",
         p_site !=1) %>% 
  mutate(pause_type = as_factor(pause_type)) %>%
  mutate(pause_type = fct_relevel(pause_type, c("induced", "resolved", "sustained"))) %>% 
  relocate(WT_pause, KO_pause, pause_type, .after = gene_sym) %>% 
  mutate(across(c(8:22), ~as.factor(.))) -> PH_AA

pauses %>% 
  filter(pause_type != "ns",
         p_site !=1) -> pauses

plot_pauses_groups <- function(df, pathway) {
  
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
    inner_join(region_lengths[,c("transcript", "CDS_len")], by = "transcript") %>% #merge with CDS lengths data
    mutate(relative_position = (p_site + 1)/(CDS_len/3)) %>%  #length normalise (CDS / 3 to convert nt to codons)
    ggplot(aes(x=relative_position, fill=pause_type)) + 
    geom_density(alpha = 0.7) +
    facet_wrap(~ gene_group, scales = "fixed", labeller = as_labeller(labeller_vec)) +
    scale_fill_manual(values = c("induced" = "red3", "sustained" = "grey40", "resolved" = "royalblue")) +
    theme_classic() +
    theme(aspect.ratio = 1) +
    labs(
      x = "% CDS",
      y = "Pause density",
      title = paste(paste(Tissue, "Npm1 KO"), str_replace_all(pathway, "_", " "), paste("Genes in pathway:", GIP), sep = "\n"),
      subtitle = paste("NES Protein", round(N_P, 2),
                       "\nadj pvalue Protein", signif(padj_P, 2),
                       "\nDetected Proteins:", DP,
                       "\nNES RPFs", round(N_R, 2),
                       "\nadj pvalue RPFs", signif(padj_R, 2),
                       "\nDetected in Riboseq:", DR)) +
    theme_bw() +
    theme(axis.title = element_text(size = 14),
          strip.text = element_text(size = 14),
          axis.text = element_text(size = 10)) -> density_plot
  
  pdf(file = paste0(plot_dir, paste(pathway, "RPF_Prot_pause_density_plot_fixedAx_allProp.pdf", sep = "_")), width = 10, height = 7)
  print(density_plot)
  dev.off()
  
  return(density_plot)
}

# On gene groups from the leading edges of gseas ---------
# for every pathway, do the plot of the pause distribution across CDSs by pause type

# fetch number of genes in each pathway - regardless of whether they are detected or not
source(paste0(parent_dir,"/R_scripts/read_mouse_GSEA_pathways.R"))

## read in leading edges HALLMARKS ---------
LE_RPFs <- read_tsv(file = file.path(parent_dir, "plots/fgsea/scatters/AK_Npm1KO_Hallmark_RPF_BC_LeadEdge.tsv")) 
LE_Protein <- read_tsv(file = file.path(parent_dir, "plots/fgsea/Proteomics/scatters/AK_Npm1KO_Hallmark_LeadEdge.tsv"))

HM <- data.frame(pathway = names(pathways.hallmark),
                 count = sapply(pathways.hallmark, length),
                 stringsAsFactors = FALSE)
HM %>% remove_rownames() -> HM # number of genes in each pathway

# as of 20250325 - this analysis done on all shared pathways and not only on the significant ones
LE_Protein %>% 
  left_join(LE_RPFs, by = join_by(pathway), suffix = c("_Prot", "_RPFs")) %>% 
  left_join(HM)-> LE

LE %>% pull(pathway) -> POI

# plot
lapply(POI, plot_pauses_groups, df = pauses)

## read in leading edges GOBP ---------
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

# plot
lapply(POI, plot_pauses_groups, df = pauses)
