# Nucleophosmin role in colorectal cancer
The scripts in this repository have been used in the analysis of Riboseq data from a murine colorectal cancer model upon depletion of Nucleophosmin (Npm1).

The results have been published in the manuscript titled "Nucleophosmin supports WNT-driven hyperproliferation and tumour initiation" (DOI here)

The basis of all analysis are the scripts in the main Riboseq branch of the [Bushell lab](https://github.com/Bushell-lab/Ribo-seq), developed by Dr. Joseph Waldron, Associate Scientist in the lab of Prof. Martin Bushell at CRUK Scotland Institute.

The pause site analysis script was written by Dr. Pauline Herviou, Postdoctoral Scientist in the lab of Prof. Martin Bushell at CRUK Scotland Institute.
The pause analysis plotting functions have then been integrated with the analyses of the genes belonging to the different Hallmark pathways leading edges.

## Raw data availability
A key aspect of the project was the integration of Riboseq and proteomics from the same samples.

The raw Riboseq / cytoplasmic RNA seq data has been deposited at the Gene Expression Omnibus under the accession [GSE249958](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE249958)

The mass spectrometry proteomics data have been deposited to the ProteomeXchange Consortium via the PRIDE partner repository with the dataset identifier PXD062969 and 10.6019/PXD062969

## Environments

The analysis is done within specific envirnments as described in the [Riboseq pipeline installation page](https://github.com/Bushell-lab/Ribo-seq/tree/main/Installation)

The exact composition of the environment used for these analysis is in the Riboseq_environment.yml and RNAseq_environment.yml files in the Environment subfolder

## Matching with proteomics data

One of the issue of MS data is that in some cases it is impossible to determine if the detected peptide describe an individual protein or a protein-group.

In order to match to Riboseq data, I made the decision to assign to all the gene_symbols in a protein group the same intensities and ttest differences calculated for the whole group (using a separate_longer_delim(gene_sym, delim = ";") line in the R scripts).

What this means is that some proteomics values are duplicated but we do not lose any Riboseq identified genes in exchange.

## Machines details

Shell scripts analyses were run on a terminal in Ubuntu 20.04.6 LTS

R Scripts were run on R 4.3.3 and 4.2.2 with no obvious differences in the results (other than slightly more warnings being triggered by using 4.2.2)
