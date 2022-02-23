## DESeq2 DGE of simulated F1 expression, generated from
## parental expression, where reads were pseudoaligned against
## Xbirchmanni gtf reference
## note that xbirch and xmal abundance files are real, observed values from
## kallisto while F1 abundance files are simulated
## see Scripts/expression_analysis/simulate_null_F1_expression.R
## updated I-2022

setwd("~/Box/Schumer_lab_resources/Project_files/Thermal_tolerance_projects/")

library(tximportData)
library(tximport)
library(GenomicFeatures)
library(readr)
library(DESeq2)

# specify tissue: "brain" or "liver"
tissue <- "liver"

## Read in sample files
dir <- "Data/kallisto_output/simulatedF1_kallisto_posttrim_xbirch-inferred-txtome"
samples <- read.table(file.path(dir, "TT_samples_simulated-F1.txt"), header = TRUE)

# make sure that non-continuous variables are cast as factors
samples$temp <- factor(samples$temp)

# subset data by tissue
samples <- samples[samples$tissue == tissue, ]
#samples <- samples[samples$temp == "22.5", ]

samples

files <- file.path(dir, paste(samples$file_basename,"kallisto",sep="_"), "abundance.tsv")
names(files) <- paste0(samples$sample)
files

## load transcript annotations
txdb <- makeTxDbFromGFF(file="Scripts/input_files/refs/xiphophorus_birchmanni_10x_12Sep2018_yDAA6_start-codons_formatted-for-txdb.gtf", format="gtf")

# create tx2gene table, match transcript name to geneid
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")

# get gene-level count information from kallisto counts
txi <- tximport(files, type = "kallisto", tx2gene = tx2gene)
names(txi)
head(txi$counts)


#################### DGE analysis ####################

# instantiate the DESeqDataSet
dds     <- DESeqDataSetFromTximport(txi,
                                    colData = samples,
                                    design = ~ replicate + species + temp + species:temp)

# Set Xmal, 22.5 as reference levels
dds$species <- relevel(dds$species, ref = "Xmal")
dds$temp <- relevel(dds$temp, ref = "22.5")

# Determine best fit
# by plotting or quantitatively with median absolute residual (smaller = better)
# local: brain=0.0230, liver=0.0938
dds <- DESeq(dds, fitType="local")
plotDispEsts(dds, main="Dispersion plot with local fit")
residual <- mcols(dds)$dispGeneEst - mcols(dds)$dispFit
absres<-abs(residual)
summary(absres)
# vs parametric: brain=0.0300, liver=0.1278
dds <- DESeq(dds, fitType="parametric")
plotDispEsts(dds, main="Dispersion plot with parametric fit")
residual <- mcols(dds)$dispGeneEst - mcols(dds)$dispFit
absres<-abs(residual)
summary(absres)

# chose local for consistency with real F1 expression analysis (since close enough)
dds <- DESeq(dds, fitType="local")
resultsNames(dds)

# save dds object
saveRDS(dds, file.path(dir, paste("TT-",tissue,"-xbirch-gtf_nullF1_dds.rds", sep="")))

# save vst object
vst <- vst(dds)
saveRDS(vst, file.path(dir,paste("TT-",tissue,"-xbirch-gtf_vst.rds", sep="")))

# get temperature-dependent comparisons per species
res.mal33cV22c <- lfcShrink(dds, coef="temp_33.5_vs_22.5", type="ashr")
res.bir33cV22c <- lfcShrink(dds, contrast=list( c("temp_33.5_vs_22.5","speciesXbirch.temp33.5") ), type="ashr")
res.F133cV22c  <- lfcShrink(dds, contrast=list( c("temp_33.5_vs_22.5","speciesmalxbirchF1.temp33.5") ), type="ashr")

# get species-dependent comparisons at 22C
res.22CbirVmal <- lfcShrink(dds, coef="species_Xbirch_vs_Xmal", type="ashr")
res.22CF1Vmal  <- lfcShrink(dds, coef="species_malxbirchF1_vs_Xmal", type="ashr")
res.22CbirVF1  <- lfcShrink(dds, contrast=list("species_Xbirch_vs_Xmal","species_malxbirchF1_vs_Xmal"), type="ashr")

# get species-dependent comparisons at 33C (need to relevel)
dds$species <- relevel(dds$species, ref = "Xmal")
dds$temp <- relevel(dds$temp, ref = "33.5")
dds <- DESeq(dds, fitType="local")

res.33CbirVmal <- lfcShrink(dds, coef="species_Xbirch_vs_Xmal", type="ashr")
res.33CF1Vmal  <- lfcShrink(dds, coef="species_malxbirchF1_vs_Xmal", type="ashr")
res.33CbirVF1  <- lfcShrink(dds, contrast=list("species_Xbirch_vs_Xmal","species_malxbirchF1_vs_Xmal"), type="ashr")


## note in the future could use nbinomWaldTest for releveling
#dds$temp <- relevel(dds$temp, ref = "33.5")
#dds <- nbinomWaldTest(dds) # rearranges design, coefs now use F1 as species reference

#################### START PRODUCE OUT TABLES ####################

# collect all results DESeq2 object variables starting with "res."
res_list <- lapply(ls(pattern = "^res\\."), get)
names(res_list)<- (ls(pattern = "^res\\."))
names(res_list)
#XXX
of_header   <- paste("nullF1_TT-",tissue,"-xbirch-gtf_",sep="")
# loop through all results objects, output to csv files
for (i in 1:length(res_list)) {

  ## Create output file with all values for current results object
  res <- res_list[[i]]
  # reorder results object by adj p-val (increasing)
  res_reorder <- res[order(res$padj), ]
  # merge results values with the normalized counts data for all samples
  res_out <- merge(as.data.frame(res_reorder), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
  # rename first column
  names(res_out)[1] <- "Gene"
  # output this to unique file
  write.csv(res_out, file = file.path(dir,paste(of_header,"DGE_lfc-shr_",names(res_list[i]), ".csv", sep = "")))

  ## Create output file with LFC and padj values from all created res. objects
  # grab only LFC and padj columns from results object (with relevant gene names as row.names)
  res_slim <- data.frame(res$log2FoldChange,res$padj)
  row.names(res_slim) <- row.names(res)
  names(res_slim) <- c( paste("LFC_",names(res_list[i]),sep=""),paste("padj_",names(res_list[i]),sep="") )

  # if you're on the first object, initialize resdata.shr
  if (i == 1) {
    resdata.shr <- res_slim
    row.names(resdata.shr) <- row.names(res_slim)
  }
  # otherwise merge new results values with master table of LFC and padj values
  else {
    resdata.shr <- merge(resdata.shr,as.data.frame(res_slim), by="row.names",sort=FALSE)
    row.names(resdata.shr) <- row.names(res_slim)
  }
}
resdata.shr
# merge results values with the normalized counts data for all samples
resdata.shr <- merge(resdata.shr, as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata.shr)[1] <- "Gene"
write.csv(resdata.shr, file=file.path(dir,paste(of_header,"DGE_lfc-shr_all.csv",sep="")))
