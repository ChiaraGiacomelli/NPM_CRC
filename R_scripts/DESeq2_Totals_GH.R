
# load libraries & set conflicts -----
library(DESeq2)
library(tidyverse)
library(vsn)
library(tximport)
library(scTenifoldNet)
library(PerformanceAnalytics)
library(corrplot)
library(Hmisc)
library(GGally)
library(conflicted)

# preferred conflicts
conflicts_prefer(purrr::reduce)
conflicts_prefer(dplyr::rename)
conflicts_prefer(dplyr::filter)

# read in common variables
source("common_variables.R")
parent_dir <- paste0(machine_dir, '/CGIACOME/AAJ_NPM/20230622_RiboSeq_APCKRAS/Ribo-seq-Ribo-seq2.0')
setwd(paste0(parent_dir, "/R_scripts"))

#create a variable for what the treatment is----
control <- "WT"
treatment <- "KO"
TMT <- "NPM KO"

#read in gene to transcript IDs map and rename and select ENSTM and ENSGM columns----
#this is used by DESeq2 and needs to be in this structure
tx2gene <- read_tsv(file = paste0(machine_dir, "/R11/bioinformatics_resources/FASTAs/mouse/GENCODE/vM27/transcript_info/gencode.vM27.pc_transcripts_gene_IDs.txt"), col_names = F)
tx2gene %>%
  dplyr::rename(GENEID = X1,
               TXNAME = X2) %>%
  select(TXNAME, GENEID) -> tx2gene

# Tissue type ----
Tissue = "APC-KRAS"

#read in the most abundant transcripts per gene csv file----
most_abundant_transcripts <- read_csv(file = file.path(parent_dir, "Analysis/most_abundant_transcripts/most_abundant_transcripts_IDs.csv"))

#import rsem data----
#set directory where rsem output is located
rsem_dir <- file.path(parent_dir, 'rsem')

#create a named vector of files (with path)
files <- file.path(rsem_dir, paste0(Total_sample_names, ".isoforms.results"))
names(files) <- Total_sample_names

#import data with txi
txi <- tximport(files, type="rsem", tx2gene=tx2gene)

#create a data frame with the condition/replicate information----
#you need to make sure this data frame is correct for your samples, the below creates one for a n=3 with EFT226 treatment.
sample_info <- data.frame(row.names = Total_sample_names,
                          condition = factor(c(rep(control, 4), rep(treatment, 4))),
                          replicate = factor(c(1:4,1:4)))

#print the data frame to visually check it has been made as expected
sample_info

#make a DESeq data set from imported data----
ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = sample_info,
                                   design = ~ condition + replicate)

#pre-filter to remove genes with less than an average of 10 counts across all samples----
keep <- rowMeans(counts(ddsTxi)) >= 10
table(keep)
ddsTxi <- ddsTxi[keep,]

#make sure levels are set appropriately so that Ctrl is "untreated"
ddsTxi$condition <- relevel(ddsTxi$condition, ref = control)

#run DESeq on DESeq data set----
dds <- DESeq(ddsTxi)

#extract results for each comparison----
res <- results(dds, contrast=c("condition", treatment, control))

#summarise results----
summary(res)

#write summary to a text file
sink(file = file.path(parent_dir, "Analysis/DESeq2_output", paste0("Totals_", Tissue, "_", treatment, "_DEseq2_summary.txt")))
summary(res)
sink()

#apply LFC shrinkage for each comparison----
lfc_shrink <- lfcShrink(dds, coef=paste("condition", treatment, "vs", control, sep = "_"), type="apeglm")

### plot MA on the spot to check for problemns -----
# if anything weird happens with apeglm shrinkage, consider other shrinking methods
pdf(file = file.path(parent_dir, "plots/DE_analysis", "MA_totals_unfiltered.pdf"))
FDR = 0.1
DESeq2::plotMA(res, alpha = FDR, main = paste0("Unshrunken MA-Plot (FDR = ", FDR, ")"), ylim = c(-5,5))
DESeq2::plotMA(lfc_shrink, alpha = FDR, main = paste0("apeglm MA-Plot (FDR = ", FDR, ")"), ylim = c(-5,5))
dev.off()

# write results to csv----
as.data.frame(lfc_shrink[order(lfc_shrink$padj),]) %>%
  rownames_to_column("gene") %>%
  inner_join(most_abundant_transcripts, by = "gene") -> DEseq2_output
write_csv(DEseq2_output, file = file.path(parent_dir, "Analysis/DESeq2_output", paste0("Totals_", Tissue, "_", treatment, "_batchCorrected_DEseq2_apeglm_LFC_shrinkage.csv")))

as.data.frame(res[order(res$padj),]) %>%
  rownames_to_column("gene") %>%
  inner_join(most_abundant_transcripts, by = "gene") -> DEseq2_output
write_csv(DEseq2_output, file = file.path(parent_dir, "Analysis/DESeq2_output", paste0("Totals_", Tissue, "_", treatment, "_batchCorrected_DEseq2_unshrunk.csv")))

# extract normalised counts and plot SD vs mean----
ntd <- normTransform(dds) #this gives log2(n + 1)
vsd <- vst(dds, blind=FALSE) #Variance stabilizing transformation
rld <- rlog(dds, blind=FALSE) #Regularized log transformation

meanSdPlot(assay(ntd))
meanSdPlot(assay(vsd))
meanSdPlot(assay(rld))

#write out normalised counts data----
#Regularized log transformation looks preferable for this data. Check for your own data and select the appropriate one
#The aim is for the range of standard deviations to be similar across the range of abundances, i.e. for the red line to be flat
as.data.frame(assay(rld)) %>%
  rownames_to_column("gene") %>%
  inner_join(most_abundant_transcripts, by = "gene") -> normalised_counts

write_csv(normalised_counts, file = file.path(parent_dir, "Analysis/DESeq2_output", paste0("Totals_", Tissue, "_", treatment, "_normalised_counts.csv")))

#plot PCA----
pcaData <- plotPCA(rld, intgroup=c("condition", "replicate"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

png(filename = file.path(parent_dir, "plots/PCAs", paste0(Tissue, "_", treatment, "_Totals_PCA.png")), width = 400, height = 350)
ggplot(pcaData, aes(PC1, PC2, color=condition, shape=replicate)) +
  geom_point(size=3) +
  geom_text(aes(label=replicate), colour = 'black',size = 6, nudge_x = 1, vjust=1)+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  theme_bw()+
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        legend.text = element_text(size = 18),
        legend.title = element_blank(),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5))+
  scale_color_manual(values = NPM_cols) +
  ggtitle(paste(Tissue, TMT, "Totals"))
dev.off()

#apply batch correct and re-plot heatmap and PCA----
mat <- assay(rld)
mat <- limma::removeBatchEffect(mat, rld$replicate)
assay(rld) <- mat

#PCA
pcaData <- plotPCA(rld, intgroup=c("condition", "replicate"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

png(filename = file.path(parent_dir, "plots/PCAs", paste0(Tissue, "_", treatment, "_Totals_batch_corrected_PCA.png")), width = 400, height = 370)
ggplot(pcaData, aes(PC1, PC2, color=condition, shape=replicate)) +
  geom_point(size=3) +
  geom_text(aes(label=replicate), colour = 'black',size = 6, nudge_x = 1, vjust=1)+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme_bw()+
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        legend.text = element_text(size = 18),
        legend.title = element_blank(),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5))+
  scale_color_manual(values = NPM_cols) +
  ggtitle(paste(Tissue, TMT, "\nTotals batch corrected"))
dev.off()

# plot correlation between replicates GGally------
normalised_counts %>%
  as_tibble() %>% 
  select(2:9) %>% 
  rename_with(~ gsub("_", " ", .x)) %>% 
  rename_with(~ gsub("Totals", "", .x)) -> nc_corr_plot

pdf(file = file.path(parent_dir, "plots/correlation", "Correlation_totals_spearman.pdf"), width = 8, height = 8)
ggpairs(nc_corr_plot,
        lower = list(continuous = wrap("points", size = 0.5)),
        diag = list(continuous = "blank"),
        upper = list(continuous = wrap("cor", method = "spearman")),
        title = "Correlation of total cytoplasmic RNA sequencing results")
dev.off()

# check for the outliers ------
normalised_counts %>% 
  as_tibble() %>% 
  select(-1, -10) %>% 
  relocate(gene_sym) %>% 
  rename_with(~gsub("_Totals", "", .x, fixed = TRUE)) -> nc_outliers_check

col_pairs <- combn(names(nc_outliers_check[,2:9]), 2, simplify = FALSE)

# For each pair, create new columns for predicted/residual/outlier, plus a PNG plot
df_extra <- map_dfc(col_pairs, function(pair) {
  
  x_col <- pair[1]
  y_col <- pair[2]
  
  pair_suffix <- paste(x_col, y_col, sep = "_")
  
  # Extract the actual vectors
  x <- nc_outliers_check[[x_col]]
  y <- nc_outliers_check[[y_col]]
  
  # Compute correlation (handle any NA or zero-variance issues)
  r_xy <- cor(x, y, use = "complete.obs", method = "spearman")
  
  # If correlation is NA (e.g., all NA or zero variance), return NA columns
  if (is.na(r_xy)) {
    message(paste("Skipping pair:", x_col, "vs", y_col, 
                  "- correlation is NA. Returning all NA for these columns."))
    return(tibble(
      !!paste0("pred_",    pair_suffix) := rep(NA_real_, length(x)),
      !!paste0("resid_",   pair_suffix) := rep(NA_real_, length(x)),
      !!paste0("outlier_", pair_suffix) := rep(NA, length(x))
    ))
  }
  
  # Slope & intercept from correlation formula
  slope <- r_xy * (sd(y, na.rm = TRUE) / sd(x, na.rm = TRUE))
  intercept <- mean(y, na.rm = TRUE) - slope * mean(x, na.rm = TRUE)
  
  # Predicted and residual
  pred  <- slope * x + intercept
  resid <- y - pred
  
  # Flag outliers if residual > N SD -> set here the SD you want as threshold
  resid_sd <- sd(resid, na.rm = TRUE)
  outlier  <- abs(resid) > (10 * resid_sd)

  # -- Plot and save as PNG --
  tmp_plot <- tibble(x = x, y = y, outlier = outlier)
  tmp_plot %>%
    arrange(outlier == TRUE) %>%
    ggplot(aes(x, y)) +
    geom_point(aes(color = outlier), size = 2) +
    geom_abline(slope = slope, intercept = intercept, color = "blue", linewidth = 1) +
    scale_color_manual(values = c("black", "red")) +
    xlab(x_col) +
    ylab(y_col) +
    labs(
      title = paste("Pairwise Spearman Comparison:\n", x_col, "vs", y_col),
      subtitle = "Flagging outliers if residuals > 10 SD", # change here the SD thershold in the label
      color = "Outlier?"
    ) -> plote

  png(filename = file.path(parent_dir, "plots/correlation", paste0(x_col, "_vs_", y_col, "_spearman_tots_10SD.png")), width = 500, height = 400)
  print(plote)
  dev.off()

  # Return a tibble with the new columns for this pair
  tibble(
    !!paste0("pred_",    pair_suffix) := pred,
    !!paste0("resid_",   pair_suffix) := resid,
    !!paste0("outlier_", pair_suffix) := outlier
  )
})

# Bind original data side-by-side with the new columns
nc_outliers_wide <- bind_cols(nc_outliers_check, df_extra)

write_csv(nc_outliers_wide, file = file.path(parent_dir, "Analysis/DESeq2_output", "outliers_check_totals_spearman_10SD.csv"))

nc_outliers_wide %>% 
  filter(if_any(starts_with("outlier_"), ~ .x == TRUE)) %>% 
  select("gene_sym", starts_with("outlier_")) -> outliers

write_csv(outliers, file = file.path(parent_dir, "Analysis/DESeq2_output", "outliers_list_totals_spearman_10SD.csv"))

