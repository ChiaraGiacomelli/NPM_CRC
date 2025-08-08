# load libraries & set conflicts -----
library(tidyverse)
library(viridis)
library(conflicted)
library(readxl)

# preferred conflicts
conflicts_prefer(purrr::reduce)
conflicts_prefer(dplyr::rename)
conflicts_prefer(dplyr::filter)

# read in common variables
source("common_variables.R")
parent_dir <- paste0(machine_dir, '/CGIACOME/AAJ_NPM/20230622_RiboSeq_APCKRAS/Ribo-seq-Ribo-seq2.0')
setwd(paste0(parent_dir, "/R_scripts"))

#create a variable for what the treatment is----
treatment <- "KO"
TMT <- "NPM KO"

# Tissue type
Tissue <- "APC-KRAS"

#themes----
mytheme <- theme_classic()+
  theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16))

# UNSHRUNK ----------
#read in DESeq2 output----
totals <- read_csv(file = file.path(parent_dir, "Analysis/DESeq2_output", paste0("Totals_", Tissue, "_", treatment, "_batchCorrected_DEseq2_unshrunk.csv")))
RPFs <- read_csv(file = file.path(parent_dir, "Analysis/DESeq2_output", paste0("RPFs_", Tissue, "_", treatment, "_batchCorrected_DEseq2_unshrunk.csv")))
TE <- read_csv((file = file.path(parent_dir, "Analysis/DESeq2_output", paste0("TE_", Tissue, "_", treatment, "_batchCorrected_DEseq2.csv"))))

## outliers ------
out_tots <- read_csv(file = file.path(parent_dir, "Analysis/DESeq2_output", "outliers_list_totals_spearman_10SD.csv"))
out_RPFs <- read_csv(file = file.path(parent_dir, "Analysis/DESeq2_output", "outliers_list_RPFs_spearman_10SD.csv"))

## s-ISR signature ---------
# from https://doi.org/10.1038/s41586-025-08794-6
# Plasticity of the mammalian integrated stress response. Nature 2025
# Supplementary table 2
# These data was used in the response to reviewers with regards to the kind of ISR which is induced upon Npm1 KO

Up_2B5_KO <- read_excel("/Users/chiara/Downloads/4th-Table2.xlsx", sheet = "Up shEif2b5 vs shCon (<0.05)")
Dn_2B5_KO <- read_excel("/Users/chiara/Downloads/4th-Table2.xlsx", sheet = "Down shEif2b5 vs shCon (<0.05)")

Up_2B5_KO %>%
  filter(FDR < 1e-2) %>%
  separate_rows(genes, sep = ":") %>%
  distinct(genes) %>%
  pull(genes) -> Up_vct

Dn_2B5_KO %>%
  filter(FDR < 1e-2) %>%
  separate_rows(genes, sep = ":") %>%
  distinct(genes) %>%
  pull(genes) -> Dn_vct

#plot volcanos----
RPFs %>%
  filter(!(is.na(padj))) %>%
  mutate(sig = factor(case_when(padj < 0.1 ~ "*",
                                padj >= 0.1 ~ "NS"))) %>%
  arrange(sig == "*") %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj), colour = sig))+
  geom_point(alpha = 0.5, size = 2.5)+
  mytheme+
  scale_color_manual(values = c(KO_color, WT_color), labels = c("*", "ns"), name = "Sign.") +
  xlab("log2FC")+
  ylab("-log10(padj)")+
  ggtitle(paste(Tissue, TMT, "\nRPFs")) -> RPFs_volcano

png(filename = file.path(parent_dir, "plots/DE_analysis", paste0("RPFs_", Tissue, "_", treatment, "_batchCorrected_volcano_unshrunk.png")), width = 400, height = 410)
print(RPFs_volcano)
dev.off()

pdf(file = file.path(parent_dir, "plots/DE_analysis", paste0("RPFs_", Tissue, "_", treatment, "_batchCorrected_volcano_unshrunk.pdf")), width = 7, height = 7.1)
print(RPFs_volcano)
dev.off()

RPFs %>%
  filter(!(is.na(padj))) %>%
  mutate(sISR = case_when(gene_sym %in% Up_vct ~ "UP",
                          gene_sym %in% Dn_vct ~"DOWN")) %>% 
  arrange(!is.na(sISR)) %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj), colour = sISR))+
  geom_point(alpha = 0.5, size = 2.5)+
  mytheme+
  scale_color_manual(values = c("darkred", "darkblue"), name = "in s-ISR") +
  xlab("log2FC")+
  ylab("-log10(padj)")+
  ggtitle(paste(Tissue, TMT, "\nRPFs")) -> RPFs_volcano

png(filename = file.path(parent_dir, "plots/DE_analysis", paste0("RPFs_", Tissue, "_", treatment, "_batchCorrected_volcano_unshrunk_sISR.png")), width = 400, height = 410)
print(RPFs_volcano)
dev.off()

totals %>%
  filter(!(is.na(padj))) %>%
  mutate(sig = factor(case_when(padj < 0.1 ~ "*",
                                padj >= 0.1 ~ "NS"))) %>%
  arrange(sig == "*") %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj), colour = sig))+
  geom_point(alpha = 0.5, size = 2.5)+
  mytheme+
  scale_color_manual(values = c(KO_color, WT_color), labels = c("*", "ns"), name = "Sign.") +
  xlab("log2FC")+
  ylab("-log10(padj)")+
  ggtitle(paste(Tissue, TMT, "\nTotals")) -> totals_volcano

png(filename = file.path(parent_dir, "plots/DE_analysis", paste0("Totals_", Tissue, "_", treatment, "_batchCorrected_volcano_unshrunk.png")), width = 400, height = 410)
print(totals_volcano)
dev.off()

pdf(file = file.path(parent_dir, "plots/DE_analysis", paste0("Totals_", Tissue, "_", treatment, "_batchCorrected_volcano_unshrunk.pdf")), width = 7, height = 7.1)
print(totals_volcano)
dev.off()

#merge RPF with totals data----
#select apdj thresholds
TE_sig_padj <- 0.1
RPF_sig_padj <- 0.1

TE_non_sig_padj <- 0.4
RPF_non_sig_padj <- 0.4

#merged data and make groups based on RPF/Total adjusted p-values or TE adjusted p-values
RPFs %>%
  select(gene, gene_sym, log2FoldChange, padj) %>%
  rename(RPFs_log2FC = log2FoldChange,
         RPFs_padj = padj) %>%
  inner_join(totals[,c("gene", "log2FoldChange", "padj", "gene_sym")], by = c("gene", "gene_sym")) %>%
  rename(totals_log2FC = log2FoldChange,
         totals_padj = padj) %>%
  inner_join(TE[,c("gene","log2FoldChange", "padj")], by = "gene") %>%
  rename(TE_log2FC = log2FoldChange,
         TE_padj = padj) %>%
  mutate(TE_group = factor(case_when(TE_padj < TE_sig_padj & TE_log2FC < 0 ~ "TE down",
                                     TE_padj < TE_sig_padj & TE_log2FC > 0 ~ "TE up",
                                     TE_padj >= TE_non_sig_padj ~ "no change",
                                     (TE_padj >= TE_sig_padj & TE_padj <= TE_non_sig_padj) | is.na(TE_padj) ~ "NS"),
                           levels = c("TE down", "no change", "TE up", "NS"), ordered = T),
         RPFs_group = factor(case_when(RPFs_padj < RPF_sig_padj & RPFs_log2FC < 0 & totals_padj >= RPF_non_sig_padj ~ "RPFs down",
                                  RPFs_padj < RPF_sig_padj & RPFs_log2FC > 0 & totals_padj >= RPF_non_sig_padj ~ "RPFs up",
                                  RPFs_padj >= RPF_non_sig_padj & totals_padj < RPF_sig_padj & totals_log2FC < 0 ~"Totals down",
                                  RPFs_padj >= RPF_non_sig_padj & totals_padj < RPF_sig_padj & totals_log2FC > 0 ~"Totals up",
                                  RPFs_padj < RPF_sig_padj & totals_padj < RPF_sig_padj & RPFs_log2FC < 0 & totals_log2FC < 0 ~ "both down",
                                  RPFs_padj < RPF_sig_padj & totals_padj < RPF_sig_padj & RPFs_log2FC > 0 & totals_log2FC > 0 ~ "both up"))) -> merged_data

summary(merged_data)

# plot TE scatters----
# based on RPFs/totals logFC
merged_data %>%
  mutate(alpha_score = case_when(is.na(RPFs_group) ~ 0.1,
                                 !(is.na(RPFs_group)) ~ 1)) %>%
  arrange(alpha_score) %>%
  ggplot(aes(x = totals_log2FC, y = RPFs_log2FC, colour = RPFs_group, alpha = alpha_score))+
  geom_abline(lty=1, colour = "grey40")+
  geom_hline(yintercept = 0, lty=1, colour = "grey40")+
  geom_hline(yintercept = 1, lty=2, colour = "grey40")+
  geom_hline(yintercept = -1, lty=2, colour = "grey40")+
  geom_vline(xintercept = 0, lty=1, colour = "grey40")+
  geom_vline(xintercept = 1, lty=2, colour = "grey40")+
  geom_vline(xintercept = -1, lty=2, colour = "grey40")+
  geom_point(size = 2)+
  scale_alpha(guide = "none")+
  scale_color_viridis_d(na.value = "grey50", name = "Regulation\ngroup", option = "magma", begin = 0.2, end = 0.8) +
  mytheme+
  xlab("Total RNA log2FC")+
  ylab("RPFs log2FC")+
  ggtitle(paste(Tissue, TMT, "TE scatter"))+
  xlim(c(-2.5,2.5))+
  ylim(c(-2.5,2.5)) -> RPF_groups_scatter_plot

png(filename = file.path(parent_dir, "plots/DE_analysis", paste0(Tissue, "_", treatment, "_batchCorrected_RPF_groups_scatter_fixedLimits_unshrunk.png")), width = 500, height = 400)
print(RPF_groups_scatter_plot)
dev.off()

pdf(file = file.path(parent_dir, "plots/DE_analysis", paste0(Tissue, "_", treatment, "_batchCorrected_RPF_groups_scatter_fixedLimits_unshrunk.pdf")), width = 7.5, height = 6)
print(RPF_groups_scatter_plot)
dev.off()

# based on split ISR gene sets
merged_data %>%
  mutate(sISR = case_when(gene_sym %in% Up_vct ~ "UP",
                          gene_sym %in% Dn_vct ~"DOWN")) %>% 
  arrange(!is.na(sISR)) %>%
  ggplot(aes(x = totals_log2FC, y = RPFs_log2FC, colour = sISR))+
  geom_abline(lty=1, colour = "grey40")+
  geom_hline(yintercept = 0, lty=1, colour = "grey40")+
  geom_hline(yintercept = 1, lty=2, colour = "grey40")+
  geom_hline(yintercept = -1, lty=2, colour = "grey40")+
  geom_vline(xintercept = 0, lty=1, colour = "grey40")+
  geom_vline(xintercept = 1, lty=2, colour = "grey40")+
  geom_vline(xintercept = -1, lty=2, colour = "grey40")+
  geom_point(size = 2, alpha = 0.5)+
  scale_color_manual(values = c("darkred", "darkblue"), name = "in s-ISR") +
  mytheme+
  xlab("Total RNA log2FC")+
  ylab("RPFs log2FC")+
  ggtitle(paste(Tissue, TMT, "TE scatter"), subtitle = "Genes in GO signatures w/ FDR < 1e-02 upon eIF2B epsilon KD")+
  xlim(c(-2.5,2.5))+
  ylim(c(-2.5,2.5)) -> RPF_groups_scatter_plot

png(filename = file.path(parent_dir, "plots/DE_analysis", paste0(Tissue, "_", treatment, "_batchCorrected_RPF_groups_scatter_fixedLimits_unshrunk_splitISR_looser.png")), width = 500, height = 400)
print(RPF_groups_scatter_plot)
dev.off()

pdf(file = file.path(parent_dir, "plots/DE_analysis", paste0(Tissue, "_", treatment, "_batchCorrected_RPF_groups_scatter_fixedLimits_unshrunk_splitISR_looser.pdf")), width = 7.5, height = 6)
print(RPF_groups_scatter_plot)
dev.off()

#write out group sizes
write.table(file = file.path(parent_dir, "plots/DE_analysis", paste0(Tissue, "_", treatment, "_batchCorrected_RPF_groups_scatter_unshrunk.txt")), summary(merged_data$RPFs_group), col.names = F, quote = F)

#write out csv
write_csv(merged_data, file.path(parent_dir, "Analysis/DESeq2_output/merged_DESeq2_AK_NPM_KO_batchCorrected_unshrunk.csv"))
