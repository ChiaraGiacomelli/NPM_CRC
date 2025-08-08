# load libraries ----
library(tidyverse)
library(ggrepel)
library(readxl)
library(conflicted)
library(fgsea)
library(data.table)

conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::select)
conflicts_prefer(purrr::reduce)
conflicts_prefer(dplyr::rename)

# read in common variables
source("common_variables.R")

# themes----
mytheme <- theme_classic()+
  theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16))

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

# FGSEA ----
## read in pathways - perso folder ----
source("read_mouse_GSEA_pathways.R")

# create a variable for what the treatment is----
control <- "WT"
treatment <- "KO"

Tissue <- "APC-KRAS"

# set the seed to ensure reproducible results
set.seed(020588)

# themes----
mytheme <- theme_classic()+
  theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        legend.position = "none")

# functions----
run_fgsea <- function(named_vector, pathway) {
  results <- fgsea(pathways = pathway,
                   stats=named_vector,
                   minSize = 10,
                   maxSize = 1000,
                   nPermSimple = 50000)
  return(results)
}

make_plot <- function(fgsea_result, padj_threshold, title) {
  plot <- ggplot(data = fgsea_result[fgsea_result$padj < padj_threshold], aes(reorder(pathway, NES), NES)) +
    geom_col(aes(fill = padj)) +
    scale_x_discrete(
      labels = function(x) {
        # Apply your transformations
        modified_labels <- gsub("_", " ", gsub("^(HALLMARK_|GOBP_|GOMF_|GOCC_|KEGG_|KEGG_MEDICUS_)", "", x))
        # Convert to title case, make lowercase, and wrap long labels
        wrapped_labels <- str_wrap(tools::toTitleCase(tolower(modified_labels)), width = 50)  # adjust width as needed
        return(wrapped_labels)
      }
    ) +
    coord_flip() +
    labs(x="Pathway", y="Normalized Enrichment Score",
         title=title) + 
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5))
  return(plot)
}

myP <- function(x) {
  p <- as.numeric(x)
  if (p < 2.2e-16) {
    p_label <- "p adj < 2.2e-16"
  } else {
    if (p < 0.001) {
      rounded_p <- formatC(p, format = "e", digits = 2)
    } else {
      rounded_p <- round(p, digits = 3)
    }
    p_label <- paste("p adj =", rounded_p)
  }
  return(p_label)
}

plot_volcano <- function(DE_data, gsea_result, gsea_set, pathway, dir) {
  
  #extract adjusted p-value and create label from it
  pval <- gsea_result$padj[gsea_result$pathway == pathway]
  plab <- myP(pval)
  NES <- gsea_result$NES[gsea_result$pathway == pathway]
  Nlab <- base::round(NES, 2)
  
  gene_names <- gsea_set[[pathway]]
  
  DE_data %>%
    mutate(group = factor(gene_sym %in% gene_names),
           alpha_score = case_when(group == T ~ 1,
                                   group == F ~ 0.1)) %>%
    mutate(pval = base::log10(pval)) %>% 
    arrange(group) %>%
    ggplot(aes(x = -(Ttest_AK_WToverKO), y = pval_AK_WToverKO, color = group, alpha = alpha_score)) +
    geom_point()+
    scale_colour_manual(values=c("grey", "red"))+
    scale_alpha(guide = "none")+
    geom_hline(yintercept = 0, lty = 2)+
    geom_vline(xintercept = 0, lty = 2)+
    mytheme+
    xlab("Ttest difference (KO-WT)")+
    ylab("p-value (-log10)")+
    ggtitle(paste("APC-KRAS NPM1 KO \n", str_replace_all(pathway, "_", " ")), subtitle = paste0(plab, "\nNES = ", Nlab)) -> scatter_plot
  
  png(filename = file.path(parent_dir, "plots/fgsea/Proteomics/scatters", dir, paste("AK_NPM1", pathway, "volcano_plot.png", sep = "_")), width = 500, height = 500)
  print(scatter_plot)
  dev.off()
}

# Run FGSEA on APC-KRAS samples ------
### make named vectors----
TMT %>%
  mutate(ttest = -(Ttest_AK_WToverKO)) %>% 
  relocate(ttest, .after = gene_sym) %>% 
  separate_longer_delim(gene_sym, delim = ";") -> TMT_longer

TMT_longer %>%
  group_by(gene_sym) %>%
  summarise(stat = mean(ttest)) %>%
  deframe() -> named_vector

### hallmark----
# carry out fgsea
fgsea_results <- run_fgsea(named_vector, pathway = pathways.hallmark)

#### write out leading edge table -----
# can be filtered for only significant pathways if desired
df<-as.data.frame(fgsea_results)
# df %>% 
#   filter(padj < 0.05) -> df
fwrite(df,file=paste0(parent_dir,"/plots/fgsea/Proteomics/scatters/","AK_Npm1KO", '_Hallmark_LeadEdge.tsv'),sep='\t',sep2=c('',' ',''))

padj <- 0.05

# plot enriched pathways
png(filename = file.path(parent_dir, "plots/fgsea/Proteomics/AK_Npm1_Hallmarks.png"), width = 500, height = 500)
make_plot(fgsea_result = fgsea_results, padj_threshold = padj, title = paste("APC-KRAS Npm1 KO Proteomics\nGSEA Hallmark gene sets"))
dev.off()

pdf(file = file.path(parent_dir, "plots/fgsea/Proteomics/AK_Npm1_Hallmarks.pdf"), width = 10, height = 10)
make_plot(fgsea_result = fgsea_results, padj_threshold = padj, title = paste("APC-KRAS Npm1 KO Proteomics\nGSEA Hallmark gene sets"))
dev.off()

#### plot overlaid volcano----
# plot volcano
sig_pathways <- fgsea_results$pathway[fgsea_results$padj < padj]
sig_pathways <- na.omit(sig_pathways)

#create directory
if (!(dir.exists(file.path(parent_dir, "plots/fgsea/Proteomics/scatters/hallmark")))) {
  dir.create(file.path(parent_dir, "plots/fgsea/Proteomics/scatters/hallmark"))
}
lapply(sig_pathways, plot_volcano, DE_data = TMT, gsea_result = fgsea_results, gsea_set = pathways.hallmark, dir = "hallmark")

#### plot enrichments as people seem to like them -----
# Everything significant
pathways = pathways.hallmark
df %>% pull(pathway) -> PTW

filtered_pathways <- pathways[names(pathways) %in% PTW]

# Example: Loop through each filtered pathway and create a plot
for(p in names(filtered_pathways)) {
  plot_title <- gsub("_", " ", p)
  
  output_file <- file.path(parent_dir, "plots", "fgsea", "Proteomics", "Enrichment_profiles", paste0("AK_NPM1KO_Prot_", p, ".pdf"))
  
  enrichment_info <- fgsea_results[fgsea_results$pathway == p, ]
  nes <- round(enrichment_info$NES, 2)
  pval <- signif(enrichment_info$padj, 2)
  
  pdf(file = output_file, width = 12, height = 8)
    plot <- plotEnrichment(filtered_pathways[[p]], named_vector) +
    labs(title = paste0(plot_title, " Proteomics Npm1 KO"),
         subtitle = paste("NES:", nes, "\nadj p-value:", pval))
    print(plot)
    dev.off()
}

### biological processes----
# carry out fgsea
fgsea_results <- run_fgsea(named_vector, pathway = pathways.bio_processes)

#### write out leading edge table -----
df<-as.data.frame(fgsea_results)
# df %>% 
#   filter(padj < 0.5) -> df
fwrite(df,file=paste0(parent_dir,"/plots/fgsea/Proteomics/scatters/","AK_Npm1KO", '_BP_LeadEdge.tsv'),sep='\t',sep2=c('',' ',''))

padj <- 0.005

# plot enriched pathways
png(filename = file.path(parent_dir, "plots/fgsea/Proteomics/AK_Npm1_BP.png"), width = 500, height = 500)
make_plot(fgsea_result = fgsea_results, padj_threshold = padj, title = paste("APC-KRAS Npm1 KO Proteomics\nGSEA Biological Processes"))
dev.off()

#### plot overlaid volcano----
# plot volcano
sig_pathways <- fgsea_results$pathway[fgsea_results$padj < padj]
sig_pathways <- na.omit(sig_pathways)

#create directory
if (!(dir.exists(file.path(parent_dir, "plots/fgsea/Proteomics/scatters/bio_processes")))) {
  dir.create(file.path(parent_dir, "plots/fgsea/Proteomics/scatters/bio_processes"))
}
lapply(sig_pathways, plot_volcano, DE_data = TMT, gsea_result = fgsea_results, gsea_set = pathways.bio_processes, dir = "bio_processes")

## plot enrichments as people seem to like them -----
# Integrated stress response is part of GOBP
pathways = pathways.bio_processes
ISR <- pathways[["GOBP_INTEGRATED_STRESS_RESPONSE_SIGNALING"]]

png(filename = file.path(parent_dir, "plots/fgsea/Proteomics/Enrichment_profiles", "AK_NPM1KO_Prot_GOBP_ISR.png"), width = 600, height = 400)
print(plotEnrichment(ISR, named_vector)+
        labs(title="ISR Proteomics Npm1 KO",
             subtitle = "GO:BP Integrated Stress Response\n"))
dev.off()

UPR <- pathways[["GOBP_RESPONSE_TO_ENDOPLASMIC_RETICULUM_STRESS"]]
lab = fgsea_results %>% 
  filter(pathway == "GOBP_RESPONSE_TO_ENDOPLASMIC_RETICULUM_STRESS") %>% 
  select(NES, padj)

png(filename = file.path(parent_dir, "plots/fgsea/Proteomics/Enrichment_profiles", "AK_NPM1KO_Prot_GOBP_UPR.png"), width = 600, height = 400)
print(plotEnrichment(UPR, named_vector)+
        labs(title="Proteomics Npm1 KO",
             subtitle = paste0("GO:BP Response to ER Stress\n",
                               "NES = ", round(lab %>% pull(NES), 2),
                               "\npadj = ", round(lab %>% pull(padj), 2))))
dev.off()

### molecular functions----
# carry out fgsea
fgsea_results <- run_fgsea(named_vector, pathway = pathways.mol_funs)

#### write out leading edge table -----
df<-as.data.frame(fgsea_results)
# df %>% 
#   filter(padj < 0.05) -> df
fwrite(df,file=paste0(parent_dir,"/plots/fgsea/Proteomics/scatters/","AK_Npm1KO", '_MF_LeadEdge.tsv'),sep='\t',sep2=c('',' ',''))

padj <- 0.01

# plot enriched pathways
png(filename = file.path(parent_dir, "plots/fgsea/Proteomics/AK_Npm1_MolecularFunctions.png"), width = 500, height = 500)
make_plot(fgsea_result = fgsea_results, padj_threshold = padj, title = paste("APC-KRAS Npm1 KO Proteomics\nGSEA Molecular Functions"))
dev.off()

#### plot overlaid volcano----
# plot volcano
sig_pathways <- fgsea_results$pathway[fgsea_results$padj < padj]
sig_pathways <- na.omit(sig_pathways)

#create directory
if (!(dir.exists(file.path(parent_dir, "plots/fgsea/Proteomics/scatters/mol_funs")))) {
  dir.create(file.path(parent_dir, "plots/fgsea/Proteomics/scatters/mol_funs"))
}
lapply(sig_pathways, plot_volcano, DE_data = TMT, gsea_result = fgsea_results, gsea_set = pathways.mol_funs, dir = "mol_funs")

### cellular component----
# carry out fgsea
fgsea_results <- run_fgsea(named_vector, pathway = pathways.cell_comp)

#### write out leading edge table -----
df<-as.data.frame(fgsea_results)
# df %>% 
#   filter(padj < 0.05) -> df
fwrite(df,file=paste0(parent_dir,"/plots/fgsea/Proteomics/scatters/","AK_Npm1KO", '_CC_LeadEdge.tsv'),sep='\t',sep2=c('',' ',''))

padj <- 0.01

# plot enriched pathways
png(filename = file.path(parent_dir, "plots/fgsea/Proteomics/AK_Npm1_CellComp.png"), width = 500, height = 500)
make_plot(fgsea_result = fgsea_results, padj_threshold = padj, title = paste("APC-KRAS Npm1 KO Proteomics\nGSEA Cellular Components"))
dev.off()

#### plot overlaid volcano----
# plot volcano
sig_pathways <- fgsea_results$pathway[fgsea_results$padj < padj]
sig_pathways <- na.omit(sig_pathways)

#create directory
if (!(dir.exists(file.path(parent_dir, "plots/fgsea/Proteomics/scatters/cell_comp")))) {
  dir.create(file.path(parent_dir, "plots/fgsea/Proteomics/scatters/cell_comp"))
}
lapply(sig_pathways, plot_volcano, DE_data = TMT, gsea_result = fgsea_results, gsea_set = pathways.cell_comp, dir = "cell_comp")

