# Differential gene expression analysis 

These scripts were used to perform differential gene expression (DGE)
analysis on RNAseq data from Payne et al 2022 thermal stress experiment 
on X. birchmanni, X. malinche, and F1 hybrids exposed to ambient and 
high temperature treatments.

All input files used in the following analysis can be found in input_files/ .

## Dependencies

R-3.6.1

The following packages are required for each anlaysis (see install.R):
* ThermTol-DESeq2, PCA ~ DESeq2, tximport, tximportData, GenomicFeatures, readr, PCAtools, ggplot2
* RNAseq-heatmap       ~ DESeq2, genefilter, pheatmap, RColorBrewer, gplots, reshape, ggplot2
* GO-analysis          ~ GOstats, GSEABase, biomaRt
* KEGG-analysis        ~ AnnotationHub, gage, pathview
* Sim-null-F1-exp      ~ truncnorm

## Contents

* **thermtol_deseq2_xbirch-annots.R, thermtol_deseq2_xmal-annots.R**

  * *Description: Full DESeq2 differential gene expression analysis, from preprocessing
                  and converting transcript counts to gene counts to calculating
                  shrunken log fold changes for expression comparisons between groups
                  (X. malinche, X. birchmanni, and F1s) at one of the two temperatures
                  (22.5C and 33.5C) and within groups between temperatures. "xbirch-annots"
                  uses the X. birchmanni transcriptome as a reference, whereas "xmal-annots" 
                  uses the X. malinche transcriptome as a reference.*

* **TT-dge_summary_info.R** 

  * *Description: Used to summarize the DGE results from DESeq2.*

* **F1exp-profile-info_xbir-ref.R**

  * *Description: Used to label F1 expression profile for each gene. Defines high and low
                  misexpression and intermediate expression relative to parent expression 
                  (high, low, intermediate).*

* **TT-RNAseq_alltissues_PCA.R**

  * *Description: Used to generate PCA plots from all tissue (brain and liver) RNAseq and from
                  single tissue RNAseq.*

* **TT-RNAseq_heatmap.R**

  * *Description: Used to generate heatmap of genes with greatest expression variance across
                  samples/groups.*

* **Payne-et-al_TT-boxplots.R**

  * *Description: Used to generate expression boxplots, as they appear in Payne et al 2021.*

* **TT-GO_analysis.R**

  * *Description: Gene Ontology enrichment analysis to identify biological pathways enriched
                  in expression response to high temperature.*

* **TT-GO_analysis_F1-misexpression.R**

  * *Description: GO enrichment analysis of genes with high or low misexpression in F1 hybrids.*

* **TT-KEGG_analysis.R**

  * *Description: KEGG enrichment analysis to identify KEGG pathways enriched in expression
                  response to high temperature.*

* **simulate_null_F1_expression.R ** 
  * *Description: Simulate F1 RNAseq transcript counts based on observed parental
                  transcript counts.* 

* **thermtol_deseq2_xbirch-annots_F1-null-simulations.R**
  * *Description: Gene expression analysis of simulated F1 expression data pseudoaligned
                  to the Xbirchmanni reference transcriptome.*

* **TT-dge_summary_info_F1-null-simulations.R**
  * *Description: Summary of brain and liver differential expression results for simulated
                  F1 gene expression.* 

## Contact

Questions and suggested improvements are greatly appreciated - please reach out to 
the repository owner at cypayne at stanford dot edu.
