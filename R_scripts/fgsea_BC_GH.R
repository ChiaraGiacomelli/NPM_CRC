#This is written for mouse data, will need to read in human pathways and edit pathway names if to be run on human data

#load libraries----
library(tidyverse)
library(fgsea)
library(Glimma)
library(data.table)

# read in common variables & other shared info ----
source("common_variables.R")
## read in pathways - perso folder ----
source("read_mouse_GSEA_pathways.R")

most_abundant_transcripts <- read_csv(file = file.path(parent_dir, "Analysis/most_abundant_transcripts/most_abundant_transcripts_IDs.csv"))

# create a variable for what the treatment is----
control <- "WT"
treatment <- "KO"

TMT <- "NPM KO"
Tissue <- "APC-KRAS"

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

collapse_pathways <- function(fgsea_results, named_vector, fgsea_pathway, padj_threshold) {
  collapsedPathways <- collapsePathways(fgsea_results[order(pval)][padj < padj_threshold], 
                                        fgsea_pathway, named_vector)
  mainPathways <- fgsea_results[pathway %in% collapsedPathways$mainPathways][
    order(-NES), pathway]
  
  return(mainPathways)
}

extract_pathways <- function(fgsea_results, named_vectors, gsea_set, padj, groups = c("RPFs", "Totals", "TE")) {
  collapsed_pathways_list <- list()
  if ("RPFs" %in% groups) {
    collapsed_pathways_list$RPFs <- collapse_pathways(fgsea_results = fgsea_results[[1]],
                                                named_vector = named_vectors[[1]],
                                                gsea_set,
                                                padj_threshold = padj)
  }
  
  if ("Totals" %in% groups) {
    collapsed_pathways_list$Totals <- collapse_pathways(fgsea_results = fgsea_results[[2]],
                                                   named_vector = named_vectors[[2]],
                                                   gsea_set,
                                                   padj_threshold = padj)
  }
  
  if ("TE" %in% groups) {
    collapsed_pathways_list$TE <- collapse_pathways(fgsea_results = fgsea_results[[3]],
                                               named_vector = named_vectors[[3]],
                                               gsea_set,
                                               padj_threshold = padj)
  }
  
  return(collapsed_pathways_list)
}

make_plot <- function(fgsea_result, padj_threshold, title) {
  # Filter the results based on padj threshold
  significant_results <- fgsea_result[fgsea_result$padj < padj_threshold, ]
  
  # If there are no significant results, create an empty plot with a message
  if (nrow(significant_results) == 0) {
    plot <- ggplot() +
      theme_void() +
      geom_text(aes(0, 0, label = "no significant results"), size = 6, hjust = 0.5) +
      labs(title = title) +
      theme(plot.title = element_text(hjust = 0.5))
    return(plot)
  }
  
  # If there are significant results, create the usual plot
  plot <- ggplot(data = significant_results, aes(reorder(pathway, NES), NES)) +
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
    labs(x = "Pathway", y = "Normalized Enrichment Score", title = title) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  
  return(plot)
}

plot_scatters <- function(df, gsea_set, pathway, dir) {
  gene_names <- gsea_set[[pathway]]
  
  df %>%
    filter(!is.na(RPFs_padj) & !is.na(totals_padj)) %>% 
    mutate(group = factor(gene_sym %in% gene_names),
           alpha_score = case_when(group == T ~ 1,
                                   group == F ~ 0.1)) %>%
    arrange(group) %>%
    ggplot(aes(x = totals_log2FC, y = RPFs_log2FC, colour = group, alpha = alpha_score))+
    geom_point()+
    scale_colour_manual(values=c("grey34", "red3"))+
    scale_alpha(guide = "none")+
    geom_abline(lty = 2)+
    geom_hline(yintercept = 0, lty = 2)+
    geom_vline(xintercept = 0, lty = 2)+
    ylim(c(-2.5,2.5))+
    xlim(c(-2.5,2.5))+
    mytheme+
    xlab("Total RNA log2FC")+
    ylab("RPFs log2FC")+
    ggtitle(paste(Tissue, TMT, str_replace_all(pathway, "_", " "), sep = "\n")) -> scatter_plot
  
  png(filename = file.path(parent_dir, "plots/fgsea/Test_publication/scatters", dir, paste(Tissue, TMT, pathway, "TE_scatter_plot_batchCorrected.png", sep = "_")), width = 500, height = 500)
  print(scatter_plot)
  dev.off()
}

## read in DESeq2 output----
DESeq2_data <- read_csv(file = file.path(parent_dir, "Analysis/DESeq2_output/merged_DESeq2_AK_NPM_KO_batchCorrected.csv"))

## make named vectors----
DESeq2_data %>%
  group_by(gene_sym) %>%
  summarise(stat = mean(RPFs_log2FC)) %>%
  deframe() -> RPFs_named_vector

DESeq2_data %>%
  group_by(gene_sym) %>%
  summarise(stat = mean(totals_log2FC)) %>%
  deframe() -> totals_named_vector

DESeq2_data %>%
  group_by(gene_sym) %>%
  summarise(stat = mean(TE_log2FC)) %>%
  deframe() -> TE_named_vector

named_vectors <- list(RPFs_named_vector, totals_named_vector, TE_named_vector)

## hallmark----
#carry out fgsea
hallmark_results <- lapply(named_vectors, run_fgsea, pathway = pathways.hallmark)

#save results
save(file = file.path(parent_dir, "Analysis/fgsea/Test_publication/hallmark_results_batchCorrected.Rdata"), hallmark_results)

#set adjusted p-value
padj <- 0.05

#plot enriched pathways
png(filename = file.path(parent_dir, "plots/fgsea/Test_publication", paste(Tissue, TMT, "RPFs_hallmark_batchCorrected.png", sep = "_")), width = 500, height = 500)
make_plot(fgsea_result = hallmark_results[[1]], padj_threshold = padj, title = paste(Tissue, TMT, "RPFs\nGSEA Hallmark gene sets"))
dev.off()

png(filename = file.path(parent_dir, "plots/fgsea/Test_publication", paste(Tissue, TMT, "totals_hallmark_batchCorrected.png", sep = "_")), width = 500, height = 500)
make_plot(fgsea_result = hallmark_results[[2]], padj_threshold = padj, title = paste(Tissue, TMT, "Total RNA\nGSEA Hallmark gene sets"))
dev.off()

png(filename = file.path(parent_dir, "plots/fgsea/Test_publication", paste(Tissue, TMT, "TE_hallmark_batchCorrected.png", sep = "_")), width = 500, height = 300)
make_plot(fgsea_result = hallmark_results[[3]], padj_threshold = padj, title = paste(Tissue, TMT, "TE\nGSEA Hallmark gene sets"))
dev.off()


pdf(file = file.path(parent_dir, "plots/fgsea/Test_publication", paste(Tissue, TMT, "RPFs_hallmark_batchCorrected.pdf", sep = "_")), width = 8, height = 6)
make_plot(fgsea_result = hallmark_results[[1]], padj_threshold = padj, title = paste(Tissue, TMT, "RPFs\nGSEA Hallmark gene sets"))
dev.off()

pdf(file = file.path(parent_dir, "plots/fgsea/Test_publication", paste(Tissue, TMT, "totals_hallmark_batchCorrected.pdf", sep = "_")), width = 8, height = 6)
make_plot(fgsea_result = hallmark_results[[2]], padj_threshold = padj, title = paste(Tissue, TMT, "Total RNA\nGSEA Hallmark gene sets"))
dev.off()

pdf(file = file.path(parent_dir, "plots/fgsea/Test_publication", paste(Tissue, TMT, "TE_hallmark_batchCorrected.pdf", sep = "_")), width = 8, height = 6)
make_plot(fgsea_result = hallmark_results[[3]], padj_threshold = padj, title = paste(Tissue, TMT, "TE\nGSEA Hallmark gene sets"))
dev.off()

#extract pathways
all_pathways <- extract_pathways(fgsea_results = hallmark_results, named_vectors = named_vectors, gsea_set = pathways.hallmark, padj = padj, groups = "RPFs")

#plot overlaid scatters
if (!(dir.exists(file.path(parent_dir, "plots/fgsea/Test_publication/scatters/hallmark")))) {
  dir.create(file.path(parent_dir, "plots/fgsea/Test_publication/scatters/hallmark"))
}
lapply(all_pathways$TE, plot_scatters, df = DESeq2_data, gsea_set = pathways.hallmark, dir = "hallmark")

#### write out leading edge table -----
df<-as.data.frame(hallmark_results[[1]])
fwrite(df,file=paste0(parent_dir,"/plots/fgsea/Test_publication/scatters/","AK_Npm1KO", '_Hallmark_RPF_BC_LeadEdge.tsv'),sep='\t',sep2=c('',' ',''))

df<-as.data.frame(hallmark_results[[2]])
fwrite(df,file=paste0(parent_dir,"/plots/fgsea/Test_publication/scatters/","AK_Npm1KO", '_Hallmark_totals_BC_LeadEdge.tsv'),sep='\t',sep2=c('',' ',''))

df<-as.data.frame(hallmark_results[[3]])
fwrite(df,file=paste0(parent_dir,"/plots/fgsea/Test_publication/scatters/","AK_Npm1KO", '_Hallmark_TE_BC_LeadEdge.tsv'),sep='\t',sep2=c('',' ',''))

#### plot enrichments as people seem to like them -----
# Everything significant at RPF level
pathways = pathways.hallmark
df<-as.data.frame(hallmark_results[[1]])
df %>%
  #filter(padj < 0.05) %>% 
  pull(pathway) -> PTW

filtered_pathways <- pathways[names(pathways) %in% PTW]
fgsea_results <- hallmark_results[[1]]

# Example: Loop through each filtered pathway and create a plot
for(p in names(filtered_pathways)) {
  plot_title <- gsub("_", " ", p)
  output_file <- file.path(parent_dir, "plots", "fgsea", "Test_publication", "Enrichment_profiles", paste0("AK_NPM1KO_RPFs_BC_", p, ".pdf"))
  
  enrichment_info <- fgsea_results[fgsea_results$pathway == p, ]
  nes <- round(enrichment_info$NES, 2)
  pval <- signif(enrichment_info$padj, 2)
  
  pdf(file = output_file, width = 12, height = 8)
  plot <- plotEnrichment(filtered_pathways[[p]], RPFs_named_vector) +
    labs(title = paste0(plot_title, " Ribosome Occupancy (RPFs) Npm1 KO"),
         subtitle = paste("NES:", nes, "\nadj p-value:", pval))
  print(plot)
  dev.off()
}

## biological processes----
#carry out fgsea
bio_processes_results <- lapply(named_vectors, run_fgsea, pathway = pathways.bio_processes)

#save results
save(file = file.path(parent_dir, "Analysis/fgsea/Test_publication/bio_processes_results_batchCorrected.Rdata"), bio_processes_results)

#### write out leading edge table -----
df<-as.data.frame(bio_processes_results[[1]])
fwrite(df,file=paste0(parent_dir,"/plots/fgsea/Test_publication/scatters/","AK_Npm1KO", '_GOBP_RPF_BC_LeadEdge.tsv'),sep='\t',sep2=c('',' ',''))

df<-as.data.frame(bio_processes_results[[2]])
fwrite(df,file=paste0(parent_dir,"/plots/fgsea/Test_publication/scatters/","AK_Npm1KO", '_GOBP_totals_BC_LeadEdge.tsv'),sep='\t',sep2=c('',' ',''))

df<-as.data.frame(bio_processes_results[[3]])
fwrite(df,file=paste0(parent_dir,"/plots/fgsea/Test_publication/scatters/","AK_Npm1KO", '_GOBP_TE_BC_LeadEdge.tsv'),sep='\t',sep2=c('',' ',''))

#set adjusted p-value
padj <- 0.01

#plot enriched pathways
png(filename = file.path(parent_dir, "plots/fgsea/Test_publication", paste(Tissue, TMT, "RPFs_bio_processes_batchCorrected.png", sep = "_")), width = 800, height = 1000)
make_plot(fgsea_result = bio_processes_results[[1]], padj_threshold = padj, title = paste(Tissue, TMT, "RPFs\nGSEA Biological Processes gene sets"))
dev.off()

png(filename = file.path(parent_dir, "plots/fgsea/Test_publication", paste(Tissue, TMT, "totals_bio_processes_batchCorrected.png", sep = "_")), width = 800, height = 1000)
make_plot(fgsea_result = bio_processes_results[[2]], padj_threshold = padj, title = paste(Tissue, TMT, "Total RNA\nGSEA Biological Processes gene sets"))
dev.off()

png(filename = file.path(parent_dir, "plots/fgsea/Test_publication", paste(Tissue, TMT, "TE_bio_processes_batchCorrected.png", sep = "_")), width = 800, height = 1000)
make_plot(fgsea_result = bio_processes_results[[3]], padj_threshold = padj, title = paste(Tissue, TMT, "TE\nGSEA Biological Processes gene sets"))
dev.off()

#extract pathways
all_pathways <- extract_pathways(fgsea_results = bio_processes_results, named_vectors = named_vectors, gsea_set = pathways.bio_processes, padj = padj, groups = "TE")

#plot overlaid scatters
if (!(dir.exists(file.path(parent_dir, "plots/fgsea/Test_publication/scatters/bio_processes")))) {
  dir.create(file.path(parent_dir, "plots/fgsea/Test_publication/scatters/bio_processes"))
}
lapply(all_pathways$TE, plot_scatters, df = DESeq2_data, gsea_set = pathways.bio_processes, dir = "bio_processes")

### plot enrichments as people seem to like them -----
# Integrated stress response is part of GOBP
pathways = pathways.bio_processes
ISR <- pathways[["GOBP_INTEGRATED_STRESS_RESPONSE_SIGNALING"]]

png(filename = file.path(parent_dir, "plots/fgsea/Test_publication/Enrichment_profiles", "AK_NPM1KO_RPFs_BC_GOBP_ISR.png"), width = 600, height = 400)
print(plotEnrichment(ISR, RPFs_named_vector)+
        labs(title="ISR RPFs Npm1 KO",
             subtitle = "GO:BP Integrated Stress Response\n"))
dev.off()

UPR <- pathways[["GOBP_RESPONSE_TO_ENDOPLASMIC_RETICULUM_STRESS"]]
lab = fgsea_results %>% 
  filter(pathway == "GOBP_RESPONSE_TO_ENDOPLASMIC_RETICULUM_STRESS") %>% 
  select(NES, padj)

png(filename = file.path(parent_dir, "plots/fgsea/Test_publication/Enrichment_profiles", "AK_NPM1KO_RPFs_BC_GOBP_UPR.png"), width = 600, height = 400)
print(plotEnrichment(UPR, RPFs_named_vector)+
        labs(title="RPFs Npm1 KO",
             subtitle = paste0("GO:BP Response to ER Stress\n",
                               "NES = ", round(lab %>% pull(NES), 2),
                               "\npadj = ", round(lab %>% pull(padj), 2))))
dev.off()

## molecular functions----
#carry out fgsea
mol_funs_results <- lapply(named_vectors, run_fgsea, pathway = pathways.mol_funs)

#save results
save(file = file.path(parent_dir, "Analysis/fgsea/Test_publication/mol_funs_results_batchCorrected.Rdata"), mol_funs_results)

#set adjusted p-value
padj <- 0.01

#plot enriched pathways
png(filename = file.path(parent_dir, "plots/fgsea/Test_publication", paste(Tissue, TMT, "RPFs_mol_funs_batchCorrected.png", sep = "_")), width = 800, height = 1000)
make_plot(fgsea_result = mol_funs_results[[1]], padj_threshold = padj, title = paste(Tissue, TMT, "RPFs\nGSEA Molecular Functions gene sets"))
dev.off()

png(filename = file.path(parent_dir, "plots/fgsea/Test_publication", paste(Tissue, TMT, "totals_mol_funs_batchCorrected.png", sep = "_")), width = 800, height = 1000)
make_plot(fgsea_result = mol_funs_results[[2]], padj_threshold = padj, title = paste(Tissue, TMT, "Total RNA\nGSEA Molecular Functions gene sets"))
dev.off()

png(filename = file.path(parent_dir, "plots/fgsea/Test_publication", paste(Tissue, TMT, "TE_mol_funs_batchCorrected.png", sep = "_")), width = 800, height = 1000)
make_plot(fgsea_result = mol_funs_results[[3]], padj_threshold = padj, title = paste(Tissue, TMT, "TE\nGSEA Molecular Functions gene sets"))
dev.off()

#extract pathways
all_pathways <- extract_pathways(fgsea_results = mol_funs_results, named_vectors = named_vectors, gsea_set = pathways.mol_funs, padj = padj, groups = "TE")

#plot overlaid scatters
if (!(dir.exists(file.path(parent_dir, "plots/fgsea/Test_publication/scatters/mol_funs")))) {
  dir.create(file.path(parent_dir, "plots/fgsea/Test_publication/scatters/mol_funs"))
}
lapply(all_pathways$TE, plot_scatters, df = DESeq2_data, gsea_set = pathways.mol_funs, dir = "mol_funs")

## cellular component----
#carry out fgsea
cell_comp_results <- lapply(named_vectors, run_fgsea, pathway = pathways.cell_comp)

#save results
save(file = file.path(parent_dir, "Analysis/fgsea/Test_publication/cell_comp_results_batchCorrected.Rdata"), cell_comp_results)

#set adjusted p-value
padj <- 0.05

#plot enriched pathways
png(filename = file.path(parent_dir, "plots/fgsea/Test_publication", paste(Tissue, TMT, "RPFs_cell_comp_batchCorrected.png", sep = "_")), width = 1000, height = 1200)
make_plot(fgsea_result = cell_comp_results[[1]], padj_threshold = padj, title = paste(Tissue, TMT, "RPFs\nGSEA Cellular Component gene sets"))
dev.off()

png(filename = file.path(parent_dir, "plots/fgsea/Test_publication", paste(Tissue, TMT, "totals_cell_comp_batchCorrected.png", sep = "_")), width = 1000, height = 1200)
make_plot(fgsea_result = cell_comp_results[[2]], padj_threshold = padj, title = paste(Tissue, TMT, "Total RNA\nGSEA Cellular Component gene sets"))
dev.off()

png(filename = file.path(parent_dir, "plots/fgsea/Test_publication", paste(Tissue, TMT, "TE_cell_comp_batchCorrected.png", sep = "_")), width = 1000, height = 800)
make_plot(fgsea_result = cell_comp_results[[3]], padj_threshold = padj, title = paste(Tissue, TMT, "TE\nGSEA Cellular Component gene sets"))
dev.off()

#extract pathways
all_pathways <- extract_pathways(fgsea_results = cell_comp_results, named_vectors = named_vectors, gsea_set = pathways.cell_comp, padj = padj, groups = "TE")

#plot overlaid scatters
if (!(dir.exists(file.path(parent_dir, "plots/fgsea/Test_publication/scatters/cell_comp")))) {
  dir.create(file.path(parent_dir, "plots/fgsea/Test_publication/scatters/cell_comp"))
}

lapply(all_pathways$TE, plot_scatters, df = DESeq2_data, gsea_set = pathways.cell_comp, dir = "cell_comp")

## tumour----
#carry out fgsea
tumour_phen_onto_results <- lapply(named_vectors, run_fgsea, pathway = pathways.tumour_phen_onto)

#save results
save(file = file.path(parent_dir, "Analysis/fgsea/Test_publication/tumour_phen_onto_results_batchCorrected.Rdata"), tumour_phen_onto_results)

#set adjusted p-value
padj <- 0.05

#plot enriched pathways
png(filename = file.path(parent_dir, "plots/fgsea/Test_publication", paste(Tissue, TMT, "RPFs_tumour_phen_onto_batchCorrected.png", sep = "_")), width = 1000, height = 700)
make_plot(fgsea_result = tumour_phen_onto_results[[1]], padj_threshold = padj, title = paste(Tissue, TMT, "RPFs\nGSEA Tumour phenotype ontology gene sets"))
dev.off()

png(filename = file.path(parent_dir, "plots/fgsea/Test_publication", paste(Tissue, TMT, "totals_tumour_phen_onto_batchCorrected.png", sep = "_")), width = 1000, height = 700)
make_plot(fgsea_result = tumour_phen_onto_results[[2]], padj_threshold = padj, title = paste(Tissue, TMT, "Total RNA\nGSEA Tumour phenotype ontology gene sets"))
dev.off()

png(filename = file.path(parent_dir, "plots/fgsea/Test_publication", paste(Tissue, TMT, "TE_tumour_phen_onto_batchCorrected.png", sep = "_")), width = 500, height = 200)
make_plot(fgsea_result = tumour_phen_onto_results[[3]], padj_threshold = padj, title = paste(Tissue, TMT, "TE\nGSEA Tumour phenotype ontology gene sets"))
dev.off()

#extract pathways
all_pathways <- extract_pathways(fgsea_results = tumour_phen_onto_results, named_vectors = named_vectors, gsea_set = pathways.tumour_phen_onto, padj = padj, groups = "TE")

#plot overlaid scatters
if (!(dir.exists(file.path(parent_dir, "plots/fgsea/Test_publication/scatters/tumour_phen_onto")))) {
  dir.create(file.path(parent_dir, "plots/fgsea/Test_publication/scatters/tumour_phen_onto"))
}
lapply(all_pathways$TE, plot_scatters, df = DESeq2_data, gsea_set = pathways.tumour_phen_onto, dir = "tumour_phen_onto")

## Curated----
#read in pathways
#carry out fgsea
curated_results <- lapply(named_vectors, run_fgsea, pathway = pathways.curated)

#save results
save(file = file.path(parent_dir, "Analysis/fgsea/Test_publication/curated_results_batchCorrected.Rdata"), curated_results)

#set adjusted p-value
padj <- 0.05

#plot enriched pathways
png(filename = file.path(parent_dir, "plots/fgsea/Test_publication", paste(Tissue, TMT, "RPFs_curated_batchCorrected.png", sep = "_")), width = 1000, height = 2000)
make_plot(fgsea_result = curated_results[[1]], padj_threshold = padj, title = paste(Tissue, TMT, "RPFs\nGSEA Curated gene sets"))
dev.off()

png(filename = file.path(parent_dir, "plots/fgsea/Test_publication", paste(Tissue, TMT, "totals_curated_batchCorrected.png", sep = "_")), width = 1000, height = 2000)
make_plot(fgsea_result = curated_results[[2]], padj_threshold = padj, title = paste(Tissue, TMT, "Total RNA\nGSEA Curated gene sets"))
dev.off()

png(filename = file.path(parent_dir, "plots/fgsea/Test_publication", paste(Tissue, TMT, "TE_curated_batchCorrected.png", sep = "_")), width = 1000, height = 1000)
make_plot(fgsea_result = curated_results[[3]], padj_threshold = padj, title = paste(Tissue, TMT, "TE\nGSEA Curated gene sets"))
dev.off()

#extract pathways
all_pathways <- extract_pathways(fgsea_results = curated_results, named_vectors = named_vectors, gsea_set = pathways.curated, padj = padj, groups = "TE")

#plot overlaid scatters
if (!(dir.exists(file.path(parent_dir, "plots/fgsea/Test_publication/scatters/curated")))) {
  dir.create(file.path(parent_dir, "plots/fgsea/Test_publication/scatters/curated"))
}

lapply(all_pathways$TE, plot_scatters, df = DESeq2_data, gsea_set = pathways.curated, dir = "curated")

## cell_type_sig----
#carry out fgsea
cell_type_sig_results <- lapply(named_vectors, run_fgsea, pathway = pathways.cell_type_sig)

#save results
save(file = file.path(parent_dir, "Analysis/fgsea/Test_publication/cell_type_sig_results_batchCorrected.Rdata"), cell_type_sig_results)

#set adjusted p-value
padj <- 0.05

#plot enriched pathways
png(filename = file.path(parent_dir, "plots/fgsea/Test_publication", paste(Tissue, TMT, "RPFs_cell_type_sig_batchCorrected.png", sep = "_")), width = 700, height = 700)
make_plot(fgsea_result = cell_type_sig_results[[1]], padj_threshold = padj, title = paste(Tissue, TMT, "RPFs\nGSEA Cell-type signature gene sets"))
dev.off()

png(filename = file.path(parent_dir, "plots/fgsea/Test_publication", paste(Tissue, TMT, "totals_cell_type_sig_batchCorrected.png", sep = "_")), width = 700, height = 700)
make_plot(fgsea_result = cell_type_sig_results[[2]], padj_threshold = padj, title = paste(Tissue, TMT, "Total RNA\nGSEA Cell-type signature gene sets"))
dev.off()

png(filename = file.path(parent_dir, "plots/fgsea/Test_publication", paste(Tissue, TMT, "TE_cell_type_sig_batchCorrected.png", sep = "_")), width = 700, height = 700)
make_plot(fgsea_result = cell_type_sig_results[[3]], padj_threshold = padj, title = paste(Tissue, TMT, "TE\nGSEA Cell-type signature gene sets"))
dev.off()

#extract pathways
all_pathways <- extract_pathways(fgsea_results = cell_type_sig_results, named_vectors = named_vectors, gsea_set = pathways.cell_type_sig, padj = padj, groups = "TE")

#plot overlaid scatters
if (!(dir.exists(file.path(parent_dir, "plots/fgsea/Test_publication/scatters/cell_type_sig")))) {
  dir.create(file.path(parent_dir, "plots/fgsea/Test_publication/scatters/cell_type_sig"))
}
lapply(all_pathways$TE, plot_scatters, df = DESeq2_data, gsea_set = pathways.cell_type_sig, dir = "cell_type_sig")

## transcription_factors----
#carry out fgsea
transcription_factors_results <- lapply(named_vectors, run_fgsea, pathway = pathways.transcription_factors)

#save results
save(file = file.path(parent_dir, "Analysis/fgsea/Test_publication/transcription_factors_results_batchCorrected.Rdata"), transcription_factors_results)

#set adjusted p-value
padj <- 0.05

#plot enriched pathways
png(filename = file.path(parent_dir, "plots/fgsea/Test_publication", paste(Tissue, TMT, "RPFs_transcription_factors_batchCorrected.png", sep = "_")), width = 700, height = 700)
make_plot(fgsea_result = transcription_factors_results[[1]], padj_threshold = padj, title = paste(Tissue, TMT, "RPFs\nGSEA Transcription factors gene sets"))
dev.off()

png(filename = file.path(parent_dir, "plots/fgsea/Test_publication", paste(Tissue, TMT, "totals_transcription_factors_batchCorrected.png", sep = "_")), width = 700, height = 700)
make_plot(fgsea_result = transcription_factors_results[[2]], padj_threshold = padj, title = paste(Tissue, TMT, "Total RNA\nGSEA Transcription factors gene sets"))
dev.off()

png(filename = file.path(parent_dir, "plots/fgsea/Test_publication", paste(Tissue, TMT, "TE_transcription_factors_batchCorrected.png", sep = "_")), width = 700, height = 700)
make_plot(fgsea_result = transcription_factors_results[[3]], padj_threshold = padj, title = paste(Tissue, TMT, "TE\nGSEA Transcription factors gene sets"))
dev.off()

#extract pathways
all_pathways <- extract_pathways(fgsea_results = transcription_factors_results, named_vectors = named_vectors, gsea_set = pathways.transcription_factors, padj = padj, groups = "TE")

#plot overlaid scatters
if (!(dir.exists(file.path(parent_dir, "plots/fgsea/Test_publication/scatters/transcription_factors")))) {
  dir.create(file.path(parent_dir, "plots/fgsea/Test_publication/scatters/transcription_factors"))
}
lapply(all_pathways$TE, plot_scatters, df = DESeq2_data, gsea_set = pathways.transcription_factors, dir = "transcription_factors")

## Custom made signatures from OS lab ----
# carry out fgsea
results <- lapply(named_vectors, run_fgsea, pathway = pathways.custom)

#save results
save(file = file.path(parent_dir, "Analysis/fgsea/Test_publication/custom_results_batchCorrected.Rdata"), results)

#set adjusted p-value
padj <- 0.5

#plot enriched pathways
png(filename = file.path(parent_dir, "plots/fgsea/Test_publication", paste(Tissue, treatment, "RPFs_Custom_batchCorrected.png", sep = "_")), width = 600, height = 400)
make_plot(fgsea_result = results[[1]], padj_threshold = padj, title = paste(Tissue, TMT, "RPFs\nGSEA Custom gene sets"))
dev.off()

png(filename = file.path(parent_dir, "plots/fgsea/Test_publication", paste(Tissue, treatment, "Totals_Custom_batchCorrected.png", sep = "_")), width = 600, height = 400)
make_plot(fgsea_result = results[[2]], padj_threshold = padj, title = paste(Tissue, TMT, "Totals\nGSEA Custom gene sets"))
dev.off()

png(filename = file.path(parent_dir, "plots/fgsea/Test_publication", paste(Tissue, treatment, "TE_Custom_batchCorrected.png", sep = "_")), width = 600, height = 400)
make_plot(fgsea_result = results[[3]], padj_threshold = padj, title = paste(Tissue, TMT, "TE\nGSEA Custom gene sets"))
dev.off()

#extract pathways
# loosen the threshold and then extract pathways so you can plot scatters also for things which are not necessarily super significant
padj <- 0.9
all_pathways <- extract_pathways(fgsea_results = results, named_vectors = named_vectors, gsea_set = pathways.custom, padj = padj, groups = "TE")

#plot overlaid scatters
if (!(dir.exists(file.path(parent_dir, "plots/fgsea/Test_publication/scatters/Custom")))) {
  dir.create(file.path(parent_dir, "plots/fgsea/Test_publication/scatters/Custom"))
}
lapply(all_pathways$TE, plot_scatters, df = DESeq2_data, gsea_set = pathways.custom, dir = "Custom")
