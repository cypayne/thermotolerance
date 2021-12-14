## For DGE heatmap visualization
# stabilize variance across the mean with regularized log
# transformation of counts
# write only the sig genes from deseq analysis

library("DESeq2")
library("genefilter")
library("pheatmap")
library("RColorBrewer")
library("gplots")
library("reshape")
library("ggplot2")

deseq_lfc <- read.csv("input_files/TT-brain-33cV22c-xbirch-gtf_DGE_lfc-shr_all.csv")
samples <- read.table("input_files/kallisto_output/TT_samples.txt", header = TRUE)
samples <- subset(samples, samples$tissue == "brain")

# load variance stabilize deseq2 normalized count data
vst_data <- readRDS("input_files/TT-brain-xbirch-gtf_vst.rds")

# convert to dataframe
vst <- assay(vst_data)
vst <- as.data.frame(vst)
vst$Gene <- rownames(vst)

# keep genes significantly differentially expressed between temperature treatments
sig_lfc_data <- subset(deseq_lfc, padj_res.bir33cV22c <= 0.1 | padj_res.F133cV22c <= 0.1 | padj_res.mal33cV22c <= 0.1)
sigGenes <- sig_lfc_data$Gene
vst <- vst[vst$Gene %in% sigGenes,]

## method 1: genefilter, pheatmap
#topVarGenes <- head(order(-rowVars(assay(vst_data))),500)
#mat <- assay(vst_data)[ topVarGenes, ]
#mat <- mat - rowMeans(mat)
# remove Gene column
vst.sub <- subset(vst, select = -Gene)

#topVarGenes <- head(order(-rowVars(as.matrix(vst.sub))),1000)
#vst.sub <- vst.sub[ topVarGenes, ]
vst.sub <- vst.sub - rowMeans(vst.sub)
df <- as.data.frame(colData(vst_data)[,c("temp","species")])
rownames(df) <- colnames(vst.sub)

# remove outliers (i.e. genes with >5 values)
#rows_to_remove <- c("g17300","g19071", "g4141", "g16388")
#vst.sub <- vst.sub[!rownames(vst.sub) %in% rows_to_remove,]

# colors
malcol_22=rgb(0/255,0/255,175/255)
hetcol_22=rgb(100/255,0/255,175/255)
bircol_22=rgb(150/255,0/255,0/255)
annotation_colors = list( species=c(Xmal=malcol_22,malxbirchF1=hetcol_22,Xbirch=bircol_22), temp=c("22.5"="blue","33.5"="red") )

# set scale
range <- max(abs(vst.sub))

pdf("TT-brain-sigGene_heatmap.pdf", w=6, h=7)
pheatmap(vst.sub, color = colorRampPalette(c("black","darkslategrey","darkcyan","white","goldenrod1","darkgoldenrod1","darkgoldenrod"))(100), breaks = seq(-range, range, length.out = 100), annotation_col=df, annotation_colors = annotation_colors, border_color=NA, show_rownames = FALSE, annotation_names_row = FALSE)
dev.off()


# with gplots (heatmap.2)
# Colors for plots below
(mycols <- brewer.pal(8, "Dark2")[1:length(unique(samples$sample))])
library(gplots)
png("qc-heatmap-samples.png", w=1000, h=1000, pointsize=20)
heatmap.2(as.matrix(mat), key=F, trace="none",
          col=colorpanel(100, "yellow", "blue"),
          ColSideColors=c("red","blue"),
          margin=c(10, 10), main="")
dev.off()


## method 2: reshape, ggplot2 
# convert to dataframe
vst <- assay(vst_data)
vst <- as.data.frame(vst)
vst$Gene <- rownames(vst)

# keep genes significantly differentially expressed between temperature treatments
sig_lfc_data <- subset(deseq_lfc, padj_res.bir33cV22c <= 0.1 | padj_res.F133cV22c <= 0.1 | padj_res.mal33cV22c <= 0.1)
sigGenes <- sig_lfc_data$Gene
vst <- vst[vst$Gene %in% sigGenes,]

# reformat vst data frame
vst <- melt(vst, id.vars=c("Gene"))

# make heatmap
heatmap <- ggplot(vst, aes(x=variable, y=Gene, fill=value)) + geom_raster() + scale_fill_viridis(trans="sqrt") + theme(axis.text.x=element_text(angle=65, hjust=1), axis.text.y=element_blank(), axis.ticks.y=element_blank())
heatmap
