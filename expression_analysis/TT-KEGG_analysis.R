## KEGG pathway enrichment analysis ##
# using DESeq2 DGE shrunken lfc and stats
# cyp III-2020

library(AnnotationHub)
library(gage)
library(pathview)

## brain
res.brain.mal33cV22c <- read.csv('input_files/TT-brain-33cV22c-w-mito_DGE_lfc-shr_res.mal33cV22c.csv', header=TRUE)
res.brain.bir33cV22c <- read.csv('input_files/TT-brain-33cV22c-w-mito_DGE_lfc-shr_res.bir33cV22c.csv', header=TRUE)
res.brain.F133cV22c <- read.csv('input_files/TT-brain-33cV22c-w-mito_DGE_lfc-shr_res.F133cV22c.csv', header=TRUE)

## liver
res.liver.mal33cV22c <- read.csv('input_files/TT-liver-33cV22c_DGE_lfc-shr_res.mal33cV22c.csv', header=TRUE)
res.liver.bir33cV22c <- read.csv('input_files/TT-liver-33cV22c_DGE_lfc-shr_res.bir33cV22c.csv', header=TRUE)
res.liver.F133cV22c <- read.csv('input_files/TT-liver-33cV22c_DGE_lfc-shr_res.F133cV22c.csv', header=TRUE)

## for each dataset, run through the rest of this script
res <- res.brain.F133cV22c
outfile_header <- "kegg.brain.F133cV22c_"

rownames(res) <- res$Gene

# subset significant genes (padj<0.1)
res.sig <- subset(res, padj < 0.1)

# Make Ensembl annotation db
# with latest version of X. maculatus annotations
ah <- AnnotationHub()
xmaDb <- query(ah, pattern = c("Xiphophorus Maculatus", "EnsDb", 99))
xmaEdb <- xmaDb[[1]]
gns <- genes(xmaEdb)

# match symbol, Entrez ID, and annotation information to results
res.sig$symbol = mapIds(xmaEdb,
                                   keys=row.names(res.sig), 
                                   column="SYMBOL",
                                   keytype="GENEID",
                                   multiVals="first")
res.sig$entrez = mapIds(xmaEdb,
                                   keys=row.names(res.sig), 
                                   column="ENTREZID",
                                   keytype="GENEID",
                                   multiVals="first")
res.sig$name =   mapIds(xmaEdb,
                                   keys=row.names(res.sig), 
                                   column="GENENAME",
                                   keytype="GENEID",
                                   multiVals="first")

# get foldchanges, match them with ENTREZ ID
foldchanges = res.sig$log2FoldChange
names(foldchanges) = res.sig$entrez
head(foldchanges)

# generate gset from X.maculatus genome ("xma")
kegg.xma.gset <- kegg.gsets(species = "xma", id.type = "kegg", check.new=TRUE)
# remove kegg.gsets metadata
kg.xma <- kegg.gsets("xma")
kegg.xma.gset <- kg.xma$kg.sets[kg.xma$sigmet.idx]

# test for changes in just one dir (all upreg or downreg)
keggres.1d <- gage(foldchanges,gsets=kegg.xma.gset, same.dir=TRUE,set.size = c(2, 500))
# test for changes in both dir (either upreg or downreg)
keggres.2d <- gage(foldchanges,gsets=kegg.xma.gset, same.dir=FALSE,set.size = c(2, 500))
head(keggres.1d$greater)  # all upreg
head(keggres.1d$less)     # all downreg
head(keggres.2d$greater)  # either dir

write.csv(keggres.1d$greater, file=paste0(outfile_header,"1d.up.csv"))
write.csv(keggres.1d$less, file=paste0(outfile_header,"1d.down.csv"))
write.csv(keggres.2d$greater, file=paste0(outfile_header,"2d.csv"))
## bir and mal brain 33v22 (1d-up): xma04141 Protein processing in endoplasmic reticulum

# get essential genes for each class
keggres.2d.esg <- esset.grp(keggres.2d$greater,foldchanges, gsets = kegg.xma.gset, same.dir = FALSE, cutoff = 0.05, test4up = T, output = F, make.plot = F)
keggres.up.esg <- esset.grp(keggres.1d$greater,foldchanges, gsets = kegg.xma.gset, same.dir = TRUE, cutoff = 0.05, test4up = T, output = F, make.plot = F)
keggres.down.esg <- esset.grp(keggres.1d$less,foldchanges, gsets = kegg.xma.gset, same.dir = TRUE, cutoff = 0.05, test4up = T, output = F, make.plot = F)
## F1 brain:
## keggres.2d.esg
##  xma04060 Cytokine-cytokine receptor interaction:
##  "102231897" "102237769" "102230371" "102216922" "102233430" "102232631" "102227458" "102235902" "102227526" "111609789"
##  xma04621 NOD-like receptor signaling pathway:
##  "102238337" "102219145" "102229465" "102226269" "102221102" "102229818" "102235142" "102223773" "102232835"
## keggres.up.esg
##  xma04080 Neuroactive ligand-receptor interaction:  "102230049" "102226973" "102237581"
##  xma00561 Glycerolipid metabolism: "102228835" "102225809"
##  xma04512 ECM-receptor interaction: "102222515" "102221825"
## keggres.down.esg
##  xma04621 NOD-like receptor signaling pathway
##  xma04142 Lysosome                       
##  xma04145 Phagosome: "102222515"
##  xma04623 Cytosolic DNA-sensing pathway
##  xma04620 Toll-like receptor signaling pathway: "102216923"
##  xma04060 Cytokine-cytokine receptor interaction: "102231897"
##  xma04514 Cell adhesion molecules (CAMs)
##  xma04622 RIG-I-like receptor signaling pathway

# check out significant genes that are in a pathway of interest
#path_genes <- res.sig[res.sig$entrez %in% c(""), ]
path_genes <- res.sig[res.sig$entrez %in% kegg.xma.gset$`xma04060 Cytokine-cytokine receptor interaction`, ]
write.csv(path_genes, paste0(outfile_header,"xma04060-Cytokine-cytokine-receptor-interaction_genes.csv"))

## Plot sig KEGG pathways 
# Get the pathways
# downreg
down_paths <- tibble::as_tibble(data.frame(id=rownames(keggres.1d$less), keggres.1d$less))
down_paths <- as.character(down_paths[1:1,]$id) # select only sig paths (pval < 0.05, qval < 0.1)

# upreg
up_paths <- tibble::as_tibble(data.frame(id=rownames(keggres.1d$greater), keggres.1d$greater))
up_paths <- as.character(up_paths[1:1,]$id) # select only sig paths (pval < 0.05, qval < 0.1)

# Get the IDs.
keggresids <- substr(down_paths, start=1, stop=8)
keggresids

# Define plotting function
plot_pathway = function(pid) pathview(gene.data=foldchanges, pathway.id=pid, species="xma", new.signature=FALSE)

# plot multiple pathways (plots saved to disk and returns a throwaway list object)
tmp = sapply(keggresids, function(pid) pathview(gene.data=foldchanges, pathway.id=pid, species="xma"))
