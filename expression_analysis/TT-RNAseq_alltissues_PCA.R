## Create PCA plot from RNAseq data, spanning tissues (brain and liver)

library(tximportData)
library(tximport)
library(GenomicFeatures)
library(readr)
library(DESeq2)
library(PCAtools)
library(ggplot2)

## Prepare count data using DESeq2 ##

# Read in sample files
dir <- "input_files/kallisto_output/"
samples <- read.table(file.path(dir, "TT_samples.txt"), header = TRUE)

# make sure that non-continuous variables are cast as factors
samples$temp <- factor(samples$temp)

## To include both liver and brain in PCA, keep both tissues
## To include only brain in PCA, subset brain
#samples <- samples[samples$tissue == "brain" | samples$tissue == "liver" ,]
samples <- samples[samples$tissue == "brain",]

## can skip and load RDS below...
files <- file.path(dir, "kallisto_posttrim_xbirch-inferred-txtome", paste(samples$file_basename,"kallisto",sep="_"), "abundance.h5")
names(files) <- paste0(samples$sample)
files

## load transcript annotations
txdb <- makeTxDbFromGFF(file="input_files/refs/xiphophorus_birchmanni_10x_12Sep2018_yDAA6_start-codons_formatted-for-txdb.gtf", format="gtf")

# create tx2gene table, match transcript name to geneid
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")

# get gene-level count information from kallisto counts
txi <- tximport(files, type = "kallisto", tx2gene = tx2gene)

# instantiate the DESeqDataSet 
dds     <- DESeqDataSetFromTximport(txi,
                                    colData = samples,
                                    design = ~ replicate + species + temp + tissue)

# Set Xmal, 22.5 as reference levels
dds$species <- relevel(dds$species, ref = "Xmal")
dds$temp <- relevel(dds$temp, ref = "22.5")

# create dds, fit local dispersion
dds <- DESeq(dds, fitType="local")

# save dds object
#saveRDS(dds, "TT-brain-liver-gill-xbirch-gtf_dds.rds")
saveRDS(dds, "input_files/TT-brain-xbirch-gtf_dds.rds")

## create PCA ## 

# read dds
#dds <- readRDS("input_files/TT-brain-liver-xbirch-gtf_dds.rds")  # for brain+liver plot
dds <- readRDS("input_files/TT-brain-xbirch-gtf_dds.rds")         # for brain plot

# variance stabilize deseq2 normalized count data
all_vst <- vst(dds)

# DESeq2's PCA functionality automatically filters out a bunch of 
# your transcripts based on low variance (biased / supervised). 
# PCAtools is unbiased / unsupervised
x <- assay(all_vst)
colnames(x)
rownames(samples) <- samples$sample
rownames(samples)
all(colnames(x) == rownames(samples))
p <- pca(x, metadata = samples, removeVar = 0.1)
pairsplot(p)

# brain PCA
biplot(p, x="PC2", y="PC3",
       lab=NULL, 
       colby = 'temp', colkey=c('22.5'='royalblue','33.5'='red3'),
       legendPosition = 'right', legendLabSize = 12, legendIconSize = 6.0,
       shape = 'species', shapekey = c('Xmal'=15, 'malxbirchF1'=17, 'Xbirch'=8),
       gridlines.major=FALSE, gridlines.minor=FALSE,
       drawConnectors = FALSE,
       title = 'TT brain PCA',
       subtitle = 'PC1 versus PC3',
       caption = '')

# brain+liver PCA
biplot(p, x="PC1", y="PC2",
       lab=NULL, 
       colby = 'tissue', colkey=c('brain'='indianred','liver'='goldenrod'),
       legendPosition = 'right', legendLabSize = 12, legendIconSize = 6.0,
       shape = 'species', shapekey = c('Xmal'=15, 'malxbirchF1'=17, 'Xbirch'=8),
       gridlines.major=FALSE, gridlines.minor=FALSE,
       drawConnectors = FALSE,
       title = 'TT brain and liver PCA',
       subtitle = 'PC1 versus PC2',
       caption = '')

# find pc that correlates most with group of interest
eigencorplot(p,
             components = getComponents(p, 1:15),
             metavars = c('species','replicate','temp'),
             col = c('darkblue', 'blue2', 'black', 'red2', 'darkred'),
             cexCorval = 0.7,
             colCorval = 'white',
             fontCorval = 2,
             posLab = 'bottomleft',
             rotLabX = 45,
             posColKey = 'top',
             cexLabColKey = 1.5,
             scale = TRUE,
             main = 'PC1-17 group correlations',
             colFrame = 'white',
             plotRsquared = FALSE)

# with pearson r^2 and p-val adjusts
eigencorplot(p,
             components = getComponents(p, 1:15),
             metavars = c('species','replicate','temp'),
             col = c('white', 'cornsilk1', 'gold', 'forestgreen', 'darkgreen'),
             cexCorval = 1.2,
             fontCorval = 2,
             posLab = 'all',
             rotLabX = 45,
             scale = TRUE,
             main = bquote(Principal ~ component ~ Pearson ~ r^2 ~ clinical ~ correlates),
             plotRsquared = TRUE,
             corFUN = 'pearson',
             corUSE = 'pairwise.complete.obs',
             signifSymbols = c('****', '***', '**', '*', ''),
             signifCutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1))


### PCA plot with reads pseudoaligned to Xmac annotated Xbirch genome ###

## Read in sample files
dir <- "input_files/kallisto_output/"
samples <- read.table(file.path(dir, "TT_samples.txt"), header = TRUE)

# make sure that non-continuous variables are cast as factors
samples$temp <- factor(samples$temp)

samples <- samples[samples$tissue == "brain",]

#files <- file.path(dir, "kall2birch", paste(samples$file_basename,"kall2birch",sep="_"), "abundance.h5")
files <- file.path(dir, "kall2birch-posttrim-w-mito", paste(samples$file_basename,"kallisto",sep="_"), "abundance.h5")
names(files) <- paste0(samples$sample)
files

## load transcript annotations
txdb <- makeTxDbFromGFF(file="input_files/refs/Xiphophorus_maculatus_LG.Xipmac4.4.2.81_unique-tx_w-mito.gtf", format="gtf")

# create tx2gene table, match transcript name to geneid
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")

# get gene-level count information from kallisto counts
txi <- tximport(files, type = "kallisto", tx2gene = tx2gene)
names(txi)
head(txi$counts)

# instantiate the DESeqDataSet 
dds     <- DESeqDataSetFromTximport(txi,
                                    colData = samples,
                                    design = ~ replicate + species + temp)

# Set Xmal, 22.5 as reference levels
dds$species <- relevel(dds$species, ref = "Xmal")
dds$temp <- relevel(dds$temp, ref = "22.5")

# fit local dispersion
dds <- DESeq(dds, fitType="local")

rld <- vst(dds, blind=FALSE)

# set group colors
#malcol=rgb(0/255,0/255,139/255)
#hetcol=rgb(65/255,105/255,225/255)
#bircol=rgb(255/255,0/255,0/255)
malcol22=rgb(0/255,0/255,175/255)
hetcol22=rgb(100/255,0/255,175/255)
bircol22=rgb(150/255,0/255,0/255)
malcol33=rgb(50/255,150/255,255/255)
hetcol33=rgb(200/255,0/255,175/255)
bircol33=rgb(255/255,0/255,0/255)

# reorder species (for legend)
rld$species <- factor(rld$species,levels=c("Xbirch","malxbirchF1","Xmal"))
# omit mito
rld2 <- rld[c(1,1:19106)]

# build PCA with DESeq2
pca <- DESeq2::plotPCA(rld2, intgroup=c("species","tissue"))
pca <- DESeq2::plotPCA(rld2, intgroup=c("species","temp"))

pca + theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
                         panel.grid.minor = element_blank(), 
                         axis.line = element_blank()) + scale_color_manual(values = c(bircol22,bircol33,hetcol22,hetcol33,malcol22,malcol33), labels=c("X.birchmanni: 22.5C", "X.birchmanni: 33.5C", "F1: 22.5C","F1: 33.5C", "X.malinche: 22.5C","X.malinche: 33.5C")) + expand_limits(x = c(-15,11), y = c(-15, 11))
pca + theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
                         panel.grid.minor = element_blank(), 
                         axis.line = element_blank()) + scale_color_manual(values = c(bircol22,bircol33,hetcol22,hetcol33,malcol22,malcol33), labels=c("X.birchmanni: 22.5C", "X.birchmanni: 33.5C", "F1: 22.5C","F1: 33.5C", "X.malinche: 22.5C","X.malinche: 33.5C")) + expand_limits(x = c(-15,11), y = c(-15, 11))


# ALTERNATIVELY: Use PCAtools
# according to PCAtools developer ~
# DESeq2's PCA functionality automatically filters out a bunch of 
# your transcripts based on low variance (biased / supervised). 
# PCAtools is unbiased / unsupervised
x <- assay(rld)
colnames(x)
rownames(samples) <- samples$sample
rownames(samples)
all(colnames(x) == rownames(samples))
p <- pca(x, metadata = samples, removeVar = 0.1)
pairsplot(p)
biplot(p)
biplot(p, x="PC3", y="PC4",
       lab = as.character(p$metadata$sample),
       colby = 'temp', colkey=c('22.5'='royalblue','33.5'='red3'),
       legendPosition = 'right', legendLabSize = 12, legendIconSize = 6.0,
       shape = 'species', shapekey = c('Xmal'=15, 'malxbirchF1'=17, 'Xbirch'=8),
       gridlines.major=FALSE, gridlines.minor=FALSE,
       drawConnectors = FALSE,
       title = 'TT brain PCA',
       subtitle = 'PC3 versus PC5',
       caption = '')

# find pc that correlates most with group of interest
eigencorplot(p,
             components = getComponents(p, 1:15),
             metavars = c('species','replicate','temp'),
             col = c('darkblue', 'blue2', 'black', 'red2', 'darkred'),
             cexCorval = 0.7,
             colCorval = 'white',
             fontCorval = 2,
             posLab = 'bottomleft',
             rotLabX = 45,
             posColKey = 'top',
             cexLabColKey = 1.5,
             scale = TRUE,
             main = 'PC1-17 group correlations',
             colFrame = 'white',
             plotRsquared = FALSE)
# with pearson r^2 and p-val adjusts
eigencorplot(p,
             components = getComponents(p, 1:15),
             metavars = c('species','replicate','temp'),
             col = c('white', 'cornsilk1', 'gold', 'forestgreen', 'darkgreen'),
             cexCorval = 1.2,
             fontCorval = 2,
             posLab = 'all',
             rotLabX = 45,
             scale = TRUE,
             main = bquote(Principal ~ component ~ Pearson ~ r^2 ~ clinical ~ correlates),
             plotRsquared = TRUE,
             corFUN = 'pearson',
             corUSE = 'pairwise.complete.obs',
             signifSymbols = c('****', '***', '**', '*', ''),
             signifCutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1))

### END PCA of all brain samples ###

