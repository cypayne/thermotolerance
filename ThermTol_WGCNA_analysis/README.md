# ThermTol WGCNA analysis

These scripts were used to create gene co-expression networks from normalized
gene counts (processed with DESeq2 [1]) using WGCNA [2], pull enriched GO pathways in
modules correlated with a trait of interest, and permute null distributions for the
expected number of hub genes and genes with transgressive F1 expression in all modules. 
The dataset that this pipeline was developed for was tissue-level RNAseq of individuals
from two parent species and their F1 hybrids, exposed to a control and a high temperature
condition. Scripts were written to be run interactively in R. 

Used to perform thermal tolerance co-expression analysis in Payne et al 2022.

All input files used in the following analysis can be found in input_files/ .

[1] Love, MI, Huber, W, Anders, S (2014) [Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8). Genome Biology, 15:550. 

[2] Langfelder P, Horvath S (2008) [WGCNA: an R package for weighted correlation network analysis](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-559). BMC Bioinformatics, 9:559. 

## Dependencies

R-3.6.1

The following packages are required for each anlaysis (see install.R):
* ThermTol_WGCNA        ~ DESeq2, WGCNA, tximportData, tximport, GenomicFeatures, readr, rhdf5
* GO_analysis_WGCNA     ~ GOstats, GSEABase, biomaRt
* WGCNA-hub-permutation ~ dplyr

## Contents

* **GO_analysis_WGCNA-brain.R, GO_analysis_WGCNA-liver.R** 
  * *Description: Gene Ontology enrichment analysis of each WGCNA module gene set in
                  the brain and liver datasets.*

* **GO_analysis_WGCNA_brain.R, GO_analysis_WGCNA-liver.R**
  * *Description: WGCNA analysis of brain and liver gene expression data to cluster
                  genes by expression profile across samples and calculate the 
                  correlation of each module's expression trend with genotype
                  and with temperature treatment.*

* **WGCNA-hub-misexp-permutation_TT-brain.R, WGCNA-hub-misexp-permutation_TT-liver.R**
  * *Description: Permute the expected number of hub genes that are misexpressed in F1s,
                 the expected number of misexpressed genes in F1s per module, and the 
                 expected number of hub genes per module.* 

## Contact

For questions and/or comments, please email me at cypayne at stanford dot edu.
