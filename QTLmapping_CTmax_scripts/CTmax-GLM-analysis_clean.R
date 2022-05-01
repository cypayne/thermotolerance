### CTmax-GLM-analysis_clean.R
#
# Script used to visualize CTmax QTL mapped using a general linear model
# Payne et al 2022
#
# Input file:
#   Output from admixture mapping =
#     numeric geno:     perform_glm_admixture_mapping_v2_gaussian_plus-multiple-cov.R
#     categorical geno: perform_glm_admixture_mapping_v2_gaussian_plus-multiple-cov.R


## load the qqman manhattan plot function
source("Scripts/QTLmapping_CTmax_scripts/adapted_qqman.R")

## Load infile
data<-read.csv(file="Scripts/input_files/final_rqtl-run_14I2021/LTREB-ctmax-oct-2020-only_ltreb-ctmax_no-xtrm-hi_1hot-site-tanks_17tanks-with-hi_categorical-geno-glm-gaussian.tsv_chr-renamed.csv",sep=",",head=TRUE)

data_trim<-{}
data_trim$CHR<-data$chr                # chromosome
data_trim$BP<-as.numeric(data$marker)  # marker position
data_trim$P<- -log(data$geno1_pval,10) # p-value
# to plot likelihood ratio instead of p-value, uncomment the following line
#data_trim$P<- data$likelihood.diff    # likelihood difference
data_trim<-as.data.frame(data_trim)
data_trim<-na.omit(data_trim)
lod_line<- -log(0.000002,10)          # calculate LOD threshold with permutation
manhattan(data_trim)
manhattan(data_trim,genomewideline=lod_line)
manhattan(data_trim,ylim=c(0,10),cex.axis=0.9)

title("LTREB ctmax genome-wide peaks: glm gaussian map ctmax~geno+sig.site.tanks")

data_trim_sub<-data_trim[data_trim$CHR=='22',]
data_trim_sub[data_trim_sub$P>lod_line,]
manhattan(data_trim_sub,genomewideline=lod_line)
title("LTREB ctmax grp22 peak: glm gaussian map")
