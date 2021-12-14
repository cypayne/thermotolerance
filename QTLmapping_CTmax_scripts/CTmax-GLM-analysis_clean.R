### CTmax-GLM-analysis_clean.R
# 
# Script used to visualize CTmax QTL mapped using a general linear model 
# in Payne et al 2021 
#
# Input file: 
#   Output from admixture mapping =
#     numeric geno:     perform_glm_admixture_mapping_v2_gaussian_plus-multiple-cov.R
#     categorical geno: perform_glm_admixture_mapping_v2_gaussian_plus-multiple-cov.R

## Load infile
data<-read.csv(file="../Data/LTREB_qtl/LTREB-ctmax-oct2020_20k-thinned_glm-gaussian_v2-mapping-infile_no-xtrm-hi_1hot-site-tanks_13tanks-without-hi_chr-renamed_at-least-140indv-marker-filter.tsv",sep="\t",head=TRUE)

data_trim<-{}
data_trim$CHR<-data$chrom             # chromosome
data_trim$BP<-as.numeric(data$marker) # marker position
data_trim$P<- -log(data$site,10)      # p-value
# to plot likelihood ratio instead of p-value, uncomment the following line
#data_trim$P<- data$likelihood.diff   # likelihood difference
data_trim<-as.data.frame(data_trim)
data_trim<-na.omit(data_trim)
lod_line<- -log(0.000002,10)          # calculate LOD threshold with permutation
manhattan(data_trim,genomewideline=lod_line)
manhattan(data_trim,ylim=c(0,10),cex.axis=0.9)

title("LTREB ctmax genome-wide peaks: glm gaussian map ctmax~geno+sig.site.tanks")

data_trim_sub<-data_trim[data_trim$CHR=='22',]
data_trim_sub[data_trim_sub$P>lod_line,]
manhattan(data_trim_sub,genomewideline=lod_line)
title("LTREB ctmax grp22 peak: glm gaussian map")
