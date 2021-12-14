### CTmax-Rqtl-analysis_clean.R
# 
# Script used to map CTmax QTL in Payne et al 2021 
#
# Input file: 
#   Rqtl format csv file, containing phenotype columns 
#   followed by marker genotype columns. Each row represents
#   one sample.
#   Genotypes are as follows: 
#     AA = homozygous X. malinche
#     BB = homozygous X. birchmanii
#     AB = heterozygote
#
#   The following filters were applied to the infile:
#     removed individudals with hi <= 0.1 and >= 0.9
#     removed individuals with het >= 0.9
#     markers thinned by 20k bp intervals

## load library
library(qtl)

setwd("~/Box/Schumer_lab_resources/Project_files/Thermal_tolerance_projects/Data/LTREB_qtl/CTmax-QTL_first-scan_data")

###### LOAD DATA ######  

## load rqtl infile (csv) = 30244 (20k bp thinned) markers, 152 individuals, 31 phenotypes
# for the 1st QTL scan
data<-read.cross("csv",dir="","LTREB-CTmax-oct2020-ONLY_20k-thinned-genos_no-xtrm-hi_22onehot-site-tanks.rqtl.csv", na.strings=c("NA"),estimate.map=FALSE)
# for the 2nd QTL scan
# data<-read.cross("csv",dir="","LTREB-CTmax-oct2020-ONLY_20k-thinned-genos_no-xtrm-hi_22onehot-site-tanks_ScyDAA6-2113-HRSCAF-2539.8811359-onehot.rqtl.csv", na.strings=c("NA"),estimate.map=FALSE)

# check that cross type, # individuals, and # markers are correct
# 30244 (20k bp thinned) markers, 152 individuals, 31 phenotypes
# 93% genoytped, 100% phenotyped
summary(data) 


###### REMOVE INDVS & MARKERS WITH EXCESS MISSING DATA ######  

## plot: visualize missing data patterns in data
par(mfrow=c(1,2), las=1)
# visualize number of markers with data per individual
plot(ntyped(data), ylab="# typed markers", main="# genotypes per individual")
# visualize number of individuals with data per marker
plot(ntyped(data, "mar"), ylab="# typed individuals", main="# genotypes per marker")

## identify names of markers to drop, based on genotypes-per-marker plot
# drop all markers with less than ~80% genotyped individuals (152*0.8 = 120)
nt.bymar <- ntyped(data, "mar")
todrop <- names(nt.bymar[nt.bymar < 120])
# drop those markers
data_sub <- drop.markers(data, todrop)

summary(data_sub) # brings total markers 30244 --> 29652, % genotyped from 93 --> 93.5%

## plot again: visualize number of markers with data per individual
#              to decide the threshold for keeping a sample
plot(ntyped(data_sub), ylab="# typed markers", main="# genotypes per individual")
plot(ntyped(data_sub, "mar"), ylab="# typed individuals", main="# genotypes per marker")

## identify individuals to drop, based on genotypes-per-individual plot
# drop all individuals with less than ~75% genotyped markers (29652*0.75 = 22300)
data_sub <- subset(data_sub, ind=(ntyped(data_sub)>22300))

summary(data_sub) # brings total individuals 152 --> 144, % genotyped from 93.5% --> 95.5%


###### REMOVE MARKERS WITH SIGNAL OF DISTORTED SEGREGATION ######

## calculate genotype frequencies and p-values for test of departure from 1:2:1 expected ratio
gt <- geno.table(data_sub) # Note: will throw a "sex ignored" warning

## look at markers that show distortion at 5% level after Bonferroni correction for multiple tests
# drop all those with greatest distortion (i.e. p-value < 1e-10)
#   xbir-mito:45 (mitochondria marker)
#   ScyDAA6-1934-HRSCAF-2318:860 to ScyDAA6-2393-HRSCAF-2888:22195903
distorted_seg <- gt[gt$P.value < 0.05/totmar(data_sub),]
todrop <-rownames(gt[gt$P.value < 1e-10,])
data_sub.ds <- drop.markers(data_sub, todrop)

summary(data_sub.ds) # brings total markers to 29652 --> 29042, final genotype proportions = AA:24.6  AB:52.3  BB:23.1


###### PLOT GENOTYPE FREQUENCIES ######

## pull out average genotype frequency across markers by individual
g <- pull.geno(data_sub.ds)
gfreq <- apply(g, 1, function(a) table(factor(a, levels=1:3)))
gfreq <- t(t(gfreq) / colSums(gfreq))

## plot avg individual genotype frequency for each genotype
par(mfrow=c(1,3), las=1)
for(i in 1:3) plot(gfreq[i,], ylab="Genotype frequency", main=c("AA", "AB", "BB")[i], ylim=c(0,1))


###### ESTIMATE RECOMBINATION ######

## estimate recombination between each pair of markers, calculate LOD for test of rec. freq. = 1/2
# Note: this operation requires high memory and time allocation, recommend running on a cluster/server, 
#       after figuring out the above parameters
data_sub.ds <- est.rf(data_sub.ds)


###### CALCULATE GENOTYPE PROBABILITIES ######

## load est.rf results after running on cluster/server
# 1st scan
data_sub.ds <- readRDS(file = "./ltreb-qtl-data_sub.ds-est.rf.rds")
# 2nd scan
#data_sub.ds <- readRDS(file = "./ltreb-qtl-data_sub.ds-est.rf_second-scan-w-qtl-covariate.rds")

## calculate conditional genotype probabilities given multipoint marker data
data_prob <- calc.genoprob(data_sub.ds)

# 1st scan
write.cross(data_prob,file="ltreb-ctmax-qtl-filtered_22-1hot-site-tank_data-prob.csv") 
# 2nd scan
#write.cross(data_prob,file="ltreb-ctmax-qtl-filtered_22-1hot-site-tank_data-prob_second-scan-w-qtl-covariate.csv")

saveRDS(data_prob,file="ltreb-ctmax-qtl-filtered_22-1hot-site-tank_data-prob.rds")


###### SELECT A MODEL ######

## determine the appropriate model to use/covariates to include

# phenotype of interest, per sample
ctmax <- as.numeric(pull.pheno(data_prob, "ctmax"))

# potential covariates, per sample
hi <- as.numeric(pull.pheno(data_prob, "hybrid_index"))    # X.malinche ancestry proportion
het <- as.numeric(pull.pheno(data_prob, "heterozygosity")) # % heterozygous genotypes across markers
site <- as.factor(pull.pheno(data_prob, "site"))           # where samples came from (STH, STL, STM)
tank <- as.factor(pull.pheno(data_prob,"tank"))            # trial tank number
site.tank <- as.factor(pull.pheno(data_prob, "site.tank")) # site-tank combination
sex <- as.factor(pull.pheno(data_prob,"sex"))              # sample sex

## select a model by calculating AIC stepwise
#   site.tank is the only recovered covariate
model_all<-lm(ctmax~hi+het+site+tank+site.tank+sex)
selectedMod <- step(model_all)

## now step through all 22 one-hot encoded site-tanks, to remove unneccessary site-tanks
STH.1 <- as.factor(pull.pheno(data_prob, "STH.1"))
STH.2 <- as.factor(pull.pheno(data_prob, "STH.2"))
STH.3 <- as.factor(pull.pheno(data_prob, "STH.3"))
STH.4 <- as.factor(pull.pheno(data_prob, "STH.4"))
STH.7 <- as.factor(pull.pheno(data_prob, "STH.7"))
STH.8 <- as.factor(pull.pheno(data_prob, "STH.8"))
STL.1 <- as.factor(pull.pheno(data_prob, "STL.1"))
STL.2 <- as.factor(pull.pheno(data_prob, "STL.2"))
STL.3 <- as.factor(pull.pheno(data_prob, "STL.3"))
STL.4 <- as.factor(pull.pheno(data_prob, "STL.4"))
STL.5 <- as.factor(pull.pheno(data_prob, "STL.5"))
STL.6 <- as.factor(pull.pheno(data_prob, "STL.6"))
STL.7 <- as.factor(pull.pheno(data_prob, "STL.7"))
STL.8 <- as.factor(pull.pheno(data_prob, "STL.8"))
STM.1 <- as.factor(pull.pheno(data_prob, "STM.1"))
STM.2 <- as.factor(pull.pheno(data_prob, "STM.2"))
STM.3 <- as.factor(pull.pheno(data_prob, "STM.3"))
STM.4 <- as.factor(pull.pheno(data_prob, "STM.4"))
STM.5 <- as.factor(pull.pheno(data_prob, "STM.5"))
STM.6 <- as.factor(pull.pheno(data_prob, "STM.6"))
STM.7 <- as.factor(pull.pheno(data_prob, "STM.7"))
STM.8 <- as.factor(pull.pheno(data_prob, "STM.8"))

## select a model with site-tanks by calculating AIC stepwise
# 1st scan
model_all<-lm(ctmax~hi+het+sex+STH.1+STH.2+STH.3+STH.4+STH.7+STH.8+
                STL.1+STL.2+STL.3+STL.4+STL.5+STL.6+STL.7+STL.8+
                STM.1+STM.2+STM.3+STM.4+STM.5+STM.6+STM.7+STM.8)
selectedMod <- step(model_all)
summary(model_all)

# pull selected covariates (plus hybrid_index)
#  ctmax ~ hi + STH.1 + STH.2 + STH.3 + STH.4 + STH.7 + STH.8 + STL.1 + 
#  STL.2 + STL.4 + STL.5 + STL.6 + STL.7 + STM.2 + STM.4 + STM.5 + 
#  STM.6 + STM.7
covariates <- pull.pheno(data_prob, c("hybrid_index","STH.1","STH.2","STH.3","STH.4","STH.7","STH.8","STL.1","STL.2","STL.4","STL.5","STL.6","STL.7","STM.2",
                                      "STM.4","STM.5","STM.6","STM.7"))

# 2nd scan
geno <- as.factor(pull.geno(data_prob, "ScyDAA6-2113-HRSCAF-2539:8811359"))
summary(lm(ctmax ~ geno + hi + STH.1 + STH.2 + STH.3 + STH.4 + STH.7 + STH.8 + STL.1 + 
             STL.2 + STL.4 + STL.5 + STL.6 + STL.7 + STM.2 + STM.4 + STM.5 + 
             STM.6 + STM.7))
# pull selected covariates
covariates <- pull.pheno(data_prob, c("hybrid_index","STH.1","STH.2","STH.3","STH.4","STH.7","STH.8","STL.1","STL.2","STL.4","STL.5","STL.6","STL.7","STM.2",
                                      "STM.4","STM.5","STM.6","STM.7","homo_birch0.0","het1.0","homo_mal2.0"))


###### RUN SINGLE-QTL SCAN ######

## perform a single-qtl scan
scanone.hk <- scanone(data_prob, method="hk",addcovar=covariates)
plot(scanone.hk,col="black")
summary(scanone.hk)
write.table(scanone.hk,file="ltreb-only-CTmax_17-1hot-site-tanks.scanone-hk_clean.tsv",sep="\t",quote=FALSE)
#write.table(scanone.hk,file="ltreb-only-CTmax_22-1hot-site-tanks.scanone-hk_second-scan-w-qtl-covariate.tsv",sep="\t",quote=FALSE)
saveRDS(scanone.hk, file="ltreb-only-CTmax_17-1hot-site-tanks.scanone-hk_clean.rds")

## run permutations to define the LOD threshold of significance
perm.hk <- scanone(data_prob, method="hk", addcovar=covariates, n.perm=1000)
saveRDS(perm.hk, "LTREB-only-ctmax_perm_17-site-tank.hk_perm.hk_clean.rds")
#saveRDS(perm.hk, "LTREB-only-ctmax_perm_site-tank.hk.no-xtreme_1hot-site-tank_perm.hk_second-scan-w-qtl-covariate.rds")

## get genome-wide LOD thresholds using genome-scan-adjusted p-values
summary(perm.hk,alpha=c(0.05, 0.1, 0.2))
## ctmax ~ hi+17sig.site.tanks
#lod
#5%  4.72
#10% 4.33
#20% 3.95
#1.5interval=4.37-1.5
cutoff_lod = 4.33 

## 2nd scan: ctmax ~ hi + sig.site.tanks + g22_qtl_genos
#cutoff_lod = 4.37 

## look at highest LOD peak and the corresponding p-value on each chromosome
summary(scanone.hk, perms=perm.hk, pvalues=TRUE)

## look at peaks that surpass the 10% LOD threshold
summary(scanone.hk, perms=perm.hk, alpha=0.1, pvalues=TRUE)
summary(scanone.hk[scanone.hk$lod > cutoff_lod,])

# 1st scan: chr22 ScyDAA6-2113-HRSCAF-2539:8811359 ScyDAA6-2113-HRSCAF-2539 8.81 4.37 0.097
# output all markers with LOD score above 10% LOD threshold
write.table(scanone.hk[scanone.hk$lod > cutoff_lod,],file="ltreb-only-CTmax_17-1hot-site-tanks.scanone-hk_peaks-above-0.1LODcutoff_clean.tsv",sep="\t",quote=FALSE)

## plot chr22 peak
pdf("ltreb-only-CTmax_17-1hot-site-tanks_chr22-qtl.pdf",8.4,8.9)
plot(scanone.hk,chr="ScyDAA6-2113-HRSCAF-2539") # chr22
abline(h=cutoff_lod, col="red", lwd=3)
title("22@8.81")
dev.off()

# 2nd scan: chr15 interacting peak
#chr  pos lod
#ScyDAA6-5984-HRSCAF-6694:5847498 ScyDAA6-5984-HRSCAF-6694 5.85 3.5
plot(scanone.hk,chr="ScyDAA6-5984-HRSCAF-6694") # chr15
abline(h=cutoff_lod, col="red", lwd=3)
title("15@5.85")


###### LOD peak support ######

## get the 1.5-LOD support interval for significant peaks
lint = cutoff_lod - 1.5   # 2.83
lodint(scanone.hk,"ScyDAA6-2113-HRSCAF-2539",drop=1.5,expandtomarkers = TRUE)
#ScyDAA6-2113-HRSCAF-2539:7006520  ScyDAA6-2113-HRSCAF-2539  7.006520 2.713022
#ScyDAA6-2113-HRSCAF-2539:8811359  ScyDAA6-2113-HRSCAF-2539  8.811359 4.366771
#ScyDAA6-2113-HRSCAF-2539:10301012 ScyDAA6-2113-HRSCAF-2539 10.301012 2.856064

## get the 2.0-LOD support interval for significant peaks
lint2 = cutoff_lod - 2   # 2.33
lodint(scanone.hk,"ScyDAA6-2113-HRSCAF-2539",drop=2,expandtomarkers = TRUE)
#ScyDAA6-2113-HRSCAF-2539:6963205  ScyDAA6-2113-HRSCAF-2539  6.963205 2.257862
#ScyDAA6-2113-HRSCAF-2539:8811359  ScyDAA6-2113-HRSCAF-2539  8.811359 4.366771
#ScyDAA6-2113-HRSCAF-2539:12213851 ScyDAA6-2113-HRSCAF-2539 12.213851 2.163837

## get 95% Bayes credible interval
bayesint(scanone.hk,"ScyDAA6-2113-HRSCAF-2539",0.95,expandtomarkers = TRUE)
#ScyDAA6-2113-HRSCAF-2539:6984383  ScyDAA6-2113-HRSCAF-2539  6.984383 2.412447
#ScyDAA6-2113-HRSCAF-2539:8811359  ScyDAA6-2113-HRSCAF-2539  8.811359 4.366771
#ScyDAA6-2113-HRSCAF-2539:12192947 ScyDAA6-2113-HRSCAF-2539 12.192947 2.479510

#out.boot <- scanoneboot(data_prob,chr="ScyDAA6-2113-HRSCAF-2539",n.boot=1000,prob=0.95)
#summary(out.boot)
#plot(out.boot)


###### MAKE EFFECT PLOTS ######

### Payne et al 2021 effect plots

##chr22 qtl effect plot
setwd("~/Box/Schumer_lab_resources/Project_files/Thermal_tolerance_projects/Data/LTREB_qtl/CTmax-QTL_first-scan_data")
data_prob.df <- read.csv("ltreb-ctmax-qtl-filtered_22-1hot-site-tank_data-prob.csv")
qtl_genos <- data_prob.df$ScyDAA6.2113.HRSCAF.2539.8811359
qtl_genos[qtl_genos == "-"] <- NA
ctmax <- data_prob.df$ctmax
qtl22_peak_data <- data.frame(ctmax=ctmax,qtl22_geno=qtl_genos)
qtl22_peak_data <- na.omit(qtl22_peak_data)

BB<-qtl22_peak_data[qtl22_peak_data$qtl22_geno=="BB",]$ctmax
MB<-qtl22_peak_data[qtl22_peak_data$qtl22_geno=="AB",]$ctmax
MM<-qtl22_peak_data[qtl22_peak_data$qtl22_geno=="AA",]$ctmax

malcol_22=rgb(0/255,0/255,175/255)
hetcol_22=rgb(100/255,0/255,175/255)
bircol_22=rgb(150/255,0/255,0/255)

## error bar function (plot 1 SD)
error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}

pdf("ltreb-only-CTmax_17-1hot-site-tanks_chr22-qtl-effect-plot-w-points_1sd.pdf",5.5,7)

plot(1:3,c(mean(MM),mean(MB),mean(BB)),col=c(malcol_22,hetcol_22,bircol_22),ylim=c(32,38),xlab="Genotype",ylab="Normalized counts",xaxt="n",pch="-",cex=3,xlim=c(0.5,3.5))
error.bar(1:3,c(mean(MM),mean(MB),mean(BB)),c(sd(MM),sd(MB),sd(BB)),col=c(malcol_22,hetcol_22,bircol_22),lwd=2)

malcol=rgb(0/255,0/255,175/255,alpha=0.3)
hetcol=rgb(100/255,0/255,175/255,alpha=0.6)
bircol=rgb(150/255,0/255,0/255,alpha=0.6)

noise<-runif(length(BB),0.2,0.35)
points(rep(3,length(BB))+noise,BB,pch=20,cex=1.8,col=bircol)

noise<-runif(length(MB),0.2,0.35)
points(rep(2,length(MB))+noise,MB,pch=20,cex=1.8,col=hetcol)

noise<-runif(length(MM),0.2,0.35)
points(rep(1,length(MM))+noise,MM,pch=20,cex=1.8,col=malcol)

mtext(c("MM","MB","BB"),at=1:3,side=1)

dev.off()

## plot the joint effects of the chr22 QTL and chr15 peak
pdf("ltreb-only-CTmax_17-1hot-site-tanks_chr22-qtl-chr15-peak-interaction-plot.pdf",5.5,7)
effectplot(data_sim, mname2="ScyDAA6-2113-HRSCAF-2539@8.81", mname1="ScyDAA6-5984-HRSCAF-6694@4.41", geno1=c("MM","MB","BB"), geno2=c("MM","MB","BB"), var.flag = "group", main = "15@4.41 x 22@8.81", xlab = "Genotype at chr22 QTL", ylab = expression('CT'['max']*' ('*~degree*C*')'), col = effect_colors, add.legend = TRUE, legend.lab = "chr15 QTL")
dev.off()



### other effect plots                            
## Load data_prob object
data_prob <- readRDS("ltreb-ctmax-qtl-filtered_22-1hot-site-tank_data-prob.rds")

## impute missing genotypes
data_sim <- sim.geno(data_prob,step=1,n.draws=16)

## plot estimated effect of QTL: phenotype vs marker genotype
malcol=rgb(0/255,0/255,175/255)
hetcol=rgb(100/255,0/255,175/255)
bircol=rgb(150/255,0/255,0/255)
effect_colors = c(malcol,hetcol,bircol)

## plot estimated effect of QTL: phenotype vs marker genotype
# 1st scan
mar_chr22 <- find.marker(data_prob, chr='ScyDAA6-2113-HRSCAF-2539', pos=8.811359)
# simple effect plot
pdf("ltreb-only-CTmax_17-1hot-site-tanks_chr22-qtl-effect-plot.pdf",8.4,8.9)
effectplot(data_sim,mname1=mar_chr22,main='22@8.81')
dev.off()
# effect plots with individuals as points
pdf("ltreb-only-CTmax_17-1hot-site-tanks_chr22-qtl-effect-plot-w-points.pdf",5.5,7)
plotPXG(data_prob,marker=mar_chr22,main='22@8.81',infer=FALSE, ylab = expression('CT'['max']*' ('*~degree*C*')'), col = effect_colors)
dev.off()

# 2nd scan
mar_chr15 <- find.marker(data_prob, chr='ScyDAA6-5984-HRSCAF-6694', pos=4.4)
plotPXG(data_prob,marker=mar_chr15,main='15@5.85')
effectplot(data_sim,mname1=mar_chr15,main='15@4.41')
# simple effect plot and effect plot with points
pdf("ltreb-only-CTmax_17-1hot-site-tanks_chr15-qtl-effect-plot_and_w-points.pdf",8.5,6)
par(mfrow=c(1,2))
effectplot(data_sim,mname1=mar_chr15,main='15@4.41', xlab = "Genotype", ylab = expression('CT'['max']*' ('*~degree*C*')'))
plotPXG(data_prob,marker=mar_chr15,main='15@4.41',infer=FALSE, ylab = expression('CT'['max']*' ('*~degree*C*')'), col = effect_colors)

dev.off()


## plot the joint effects of the chr22 QTL and chr15 peak
pdf("ltreb-only-CTmax_17-1hot-site-tanks_chr22-qtl-chr15-peak-interaction-plot.pdf",5.5,7)
effectplot(data_sim, mname2="ScyDAA6-2113-HRSCAF-2539@8.81", mname1="ScyDAA6-5984-HRSCAF-6694@4.41", geno1=c("MM","MB","BB"), geno2=c("MM","MB","BB"), var.flag = "group", main = "15@4.41 x 22@8.81", xlab = "Genotype at chr22 QTL", ylab = expression('CT'['max']*' ('*~degree*C*')'), col = effect_colors, add.legend = TRUE, legend.lab = "chr15 QTL")
dev.off()

pdf("ltreb-only-CTmax_17-1hot-site-tanks_chr22-qtl-chr15-peak-interaction-plot-w-points.pdf",9.8,8.2)
plotPXG(data_prob, marker=c(mar_chr22, mar_chr15))
dev.off()

# to get vals
# g15 interactor
# additive: 0.2425
# dominance: -0.2835
eff <- effectplot(data_sim,mname1=mar_chr22)
eff
# ScyDAA6-2113-HRSCAF-2539:8811359
#      mean      SE
#  AA: 35.98525, 0.1777719
#  AB: 35.65896, 0.1180658
#  BB: 36.00340, 0.1713042

## plot estimated QTL effects along chr22 (shows additive and dominance effects)
# effect summary for 22@8.81: 
# a           d             se.a        se.d
# 0.008494465	-0.337265149	0.124726314	0.172637572
effect_table<-effectscan(data_sim, pheno.col="ctmax", c("ScyDAA6-2113-HRSCAF-2539"), draw=F, get.se=T)
write.table(effect_table,file="ltreb-only-CTmax_17-1hot-site-tanks.chr22-qtl_effect-table.tsv",sep="\t",quote=FALSE)
# plot effect across chr22
pdf("ltreb-only-CTmax_17-1hot-site-tanks_chr22-additive-dominance-effects-plot.pdf",9.2,6.8)
effectscan(data_sim, pheno.col="ctmax", c("ScyDAA6-2113-HRSCAF-2539"), draw=T, get.se=T)
dev.off()


###### MULTI-QTL ANALYSIS ######

# create QTL object with both loci
# what="prob" will out the QTL genotype probabilities for use in HK regression
qtl <- makeqtl(data_prob, chr=c("ScyDAA6-2113-HRSCAF-2539","ScyDAA6-5984-HRSCAF-6694"), pos=c(8.81, 4.41), what="prob")

# set covariates
covariates <- pull.pheno(data_prob, c("hybrid_index","STH.1","STH.2","STH.3","STH.4","STH.7","STH.8","STL.1","STL.2","STL.4","STL.5","STL.6","STL.7","STM.2",
                                      "STM.4","STM.5","STM.6","STM.7"))

# fit the two locus additive model
# “drop one term at a time” table compares the fit of the two-QTL model to the reduced model where one QTL is omitted.
out.fq <- fitqtl(data_prob, pheno.col=1, qtl=qtl, method="hk", get.ests=TRUE, covar=covariates, formula=y~Q1+Q2 + hybrid_index + STH.1 + STH.2 + STH.3 + STH.4 + STH.7 + STH.8 + STL.1 + STL.2 + STL.4 + STL.5 + STL.6 + STL.7 + STM.2 + STM.4 + STM.5 + STM.6 + STM.7)
summary(out.fq)

# determine whether there is an interaction between the two QTL by fitting the model with the interaction
out.fqi <- fitqtl(data_prob, pheno.col=1, qtl=qtl, method="hk", get.ests=TRUE, covar=covariates, formula=y~Q1+Q2+Q1:Q2 + hybrid_index + STH.1 + STH.2 + STH.3 + STH.4 + STH.7 + STH.8 + STL.1 + STL.2 + STL.4 + STL.5 + STL.6 + STL.7 + STM.2 + STM.4 + STM.5 + STM.6 + STM.7)
summary(out.fqi)

# also assess interaction with addint, which adds one interaction at a time
addint(data_prob, qtl=qtl, method="hk")

# refine the qtl positions
# each QTL is moved to the position giving the highest likelihood, 
# and the entire process is repeated until no further improvement in likelihood can be obtained
rqtl <- refineqtl(data_prob, qtl=qtl, method="hk")
rqtl
#                             name                      chr     pos n.gen
# Q1 ScyDAA6-2113-HRSCAF-2539@10.0 ScyDAA6-2113-HRSCAF-2539 10.0490     3
# Q2  ScyDAA6-5984-HRSCAF-6694@5.7 ScyDAA6-5984-HRSCAF-6694  5.7399     3

# look for additional qtl, using the refined qtl positions
out.aq_1qtl <- addqtl(data_prob, qtl=rqtl, method="hk", covar=covariates, formula=y~Q1 + hybrid_index + STH.1 + STH.2 + STH.3 + STH.4 + STH.7 + STH.8 + STL.1 + STL.2 + STL.4 + STL.5 + STL.6 + STL.7 + STM.2 + STM.4 + STM.5 + STM.6 + STM.7)
out.aq_2qtl <- addqtl(data_prob, qtl=rqtl, method="hk", covar=covariates, formula=y~Q1+Q2+Q1:Q2 + hybrid_index + STH.1 + STH.2 + STH.3 + STH.4 + STH.7 + STH.8 + STL.1 + STL.2 + STL.4 + STL.5 + STL.6 + STL.7 + STM.2 + STM.4 + STM.5 + STM.6 + STM.7)

plot(out.aq_1qtl)
summary(out.aq_1qtl)
out.aq_1qtl[out.aq_1qtl$lod > 3,]
# as expected, recovers chr15 qtl
# ScyDAA6-5984-HRSCAF-6694:5847498  ScyDAA6-5984-HRSCAF-6694 5.85e+00 3.503
plot(out.aq_2qtl)
summary(out.aq_2qtl)
out.aq_2qtl[out.aq_2qtl$lod > 3,]
#                                                        chr      pos   lod
# ScyDAA6-1854-HRSCAF-2213:18020066 ScyDAA6-1854-HRSCAF-2213 1.80e+01 3.512
# ScyDAA6-2469-HRSCAF-2980:1032301  ScyDAA6-2469-HRSCAF-2980 1.03e+00 3.360

# grab the 1.5-LOD for both of these
lod_cutoff = 3.512-1.5
write.csv(out.aq_2qtl[out.aq_2qtl$chr == "ScyDAA6-1854-HRSCAF-2213" & out.aq_2qtl$lod > lod_cutoff,], "CTmax-chr22-chr15-add-qtl_ScyDAA6-1854-HRSCAF-2213_1.5LOD.csv")
lod_cutoff = 3.36-1.5
write.csv(out.aq_2qtl[out.aq_2qtl$chr == "ScyDAA6-2469-HRSCAF-2980" & out.aq_2qtl$lod > lod_cutoff,], "CTmax-chr22-chr15-add-qtl_ScyDAA6-1854-HRSCAF-2213_1.5LOD.csv")

###### ROUGHLY ESTIMATE EFFECT OF QTL WITH AN LM ######

### get a very rough estimate of the QTL effect size with R^2
## load infile as data.frame and do correlations/lm with that
setwd("~/Box/Schumer_lab_resources/Project_files/Thermal_tolerance_projects/Data/LTREB_qtl/CTmax-QTL_first-scan_data")
data_prob.df <- read.csv("ltreb-ctmax-qtl-filtered_22-1hot-site-tank_data-prob.csv")

## grab all relevant variables
ctmax <- data_prob.df$ctmax
hi <- data_prob.df$hybrid_index
STH.1<-data_prob.df$STH.1
STH.2<-data_prob.df$STH.2
STH.3<-data_prob.df$STH.3
STH.4<-data_prob.df$STH.4
STH.7<-data_prob.df$STH.7
STH.8<-data_prob.df$STH.8
STL.1<-data_prob.df$STL.1
STL.2<-data_prob.df$STL.2
STL.4<-data_prob.df$STL.4
STL.5<-data_prob.df$STL.5
STL.6<-data_prob.df$STL.6
STL.7<-data_prob.df$STL.7
STM.2<-data_prob.df$STM.2
STM.4<-data_prob.df$STM.4
STM.5<-data_prob.df$STM.5
STM.6<-data_prob.df$STM.6
STM.7<-data_prob.df$STM.7

# grab the genotypes of the locus at the QTL peak
chr22_geno <- data_prob.df$ScyDAA6.2113.HRSCAF.2539.8811359
chr15_geno <- data_prob.df$ScyDAA6.5984.HRSCAF.6694.4408039

## optional: output genotypes for markers that fall under the QTL / within the 2-LOD interval: 6.963205 - 12.213851
firstcol = which(colnames(data_prob.df)=="ScyDAA6.2113.HRSCAF.2539.6963205")
lastcol = which(colnames(data_prob.df)=="ScyDAA6.2113.HRSCAF.2539.12213851")
qtl_genos <- data_prob.df[c(firstcol:lastcol)]
qtl_geno_outdf <- data.frame(ctmax,hi,STH.1,STH.2,STH.3,STH.4,STH.7,STH.8,STL.1,STL.2,STL.4,STL.5,STL.6,STL.7,STM.2,STM.4,STM.5,STM.6,STM.7,qtl_genos)
write.csv(qtl_geno_outdf,'./ltreb-ctmax-qtl_genos-all-markers-under-2LOD-chr22-qtl.csv',row.names=F,quote=F)

### chr22 QTL LM effect size estimate and p-values
## organize variables of interest into new data frame, remove missing data
dp_noNA<-data.frame(ctmax,hi,chr22_geno,STH.1,STH.2,STH.3,STH.4,STH.7,STH.8,STL.1,STL.2,STL.4,STL.5,STL.6,STL.7,STM.2,STM.4,STM.5,STM.6,STM.7)
dp_noNA[dp_noNA == '-']<-NA
dp_noNA<-na.omit(dp_noNA)
dim(dp_noNA)

## perform a linear regression with the full model and get the R^2 value 
# full model adjusted R^2: 0.6731 (pval < 2e-16)
lm.res<-lm(ctmax ~ chr22_geno + hi + STH.1 + STH.2 + STH.3 + STH.4 + STH.7 + STH.8 + STL.1 + 
             STL.2 + STL.4 + STL.5 + STL.6 + STL.7 + STM.2 + STM.4 + STM.5 + 
             STM.6 + STM.7,data=dp_noNA)
summary(lm.res)

## get the adjusted R^2 from the null model (i.e. without genotype as a coefficient)
# null adj R^2: 0.6289 (pval < 2e-16)
null.lm.res<-lm(ctmax ~ hi + STH.1 + STH.2 + STH.3 + STH.4 + STH.7 + STH.8 + STL.1 + 
                  STL.2 + STL.4 + STL.5 + STL.6 + STL.7 + STM.2 + STM.4 + STM.5 + 
                  STM.6 + STM.7,data=dp_noNA)
summary(null.lm.res)

## subtract the adjusted R^2 of the null model from the full model to get effect estimate
effect = 0.6731 - 0.6289 # 0.0442 (4.42%)


## get log likelihood difference between focal and null models
summary(lm.res)
model2<-glm(as.numeric(ctmax)~as.numeric(hi)+
              as.factor(STH.1) + as.factor(STH.2) + as.factor(STH.3) + 
              as.factor(STH.4) + as.factor(STH.7) + as.factor(STH.8) + 
              as.factor(STL.1) + as.factor(STL.2) + as.factor(STL.4) + 
              as.factor(STL.5) + as.factor(STL.6) + as.factor(STL.7) + 
              as.factor(STM.2) + as.factor(STM.4) + as.factor(STM.5) + 
              as.factor(STM.6) + as.factor(STM.7),data=dp_noNA,family="gaussian")
null<-logLik(model2)[1]

model1<-glm(as.numeric(ctmax)~as.numeric(hi)+
              as.factor(STH.1) + as.factor(STH.2) + as.factor(STH.3) + 
              as.factor(STH.4) + as.factor(STH.7) + as.factor(STH.8) + 
              as.factor(STL.1) + as.factor(STL.2) + as.factor(STL.4) + 
              as.factor(STL.5) + as.factor(STL.6) + as.factor(STL.7) + 
              as.factor(STM.2) + as.factor(STM.4) + as.factor(STM.5) + 
              as.factor(STM.6) + as.factor(STM.7) + as.factor(geno),data=dp_noNA,family="gaussian")
focal<-logLik(model1)[1]

# calculate log likelihood difference
like_diff<-focal-null
like_diff # 10.03532 


### chr15 QTL LM effect size estimate
## organize variables of interest into new data frame, remove missing data
dp_noNA<-data.frame(ctmax,hi,chr15_geno,STH.1,STH.2,STH.3,STH.4,STH.7,STH.8,STL.1,STL.2,STL.4,STL.5,STL.6,STL.7,STM.2,STM.4,STM.5,STM.6,STM.7)
dp_noNA[dp_noNA == '-']<-NA
dp_noNA<-na.omit(dp_noNA)
dim(dp_noNA)
# full model adjusted R^2: 0.6417 (pval < 2e-16)
lm.res<-lm(ctmax ~ chr15_geno + hi + STH.1 + STH.2 + STH.3 + STH.4 + STH.7 + STH.8 + STL.1 + 
             STL.2 + STL.4 + STL.5 + STL.6 + STL.7 + STM.2 + STM.4 + STM.5 + 
             STM.6 + STM.7,data=dp_noNA)
summary(lm.res)

## get the adjusted R^2 from the null model (i.e. without genotype as a coefficient)
# null adj R^2: 0.6203 (pval < 2e-16)
null.lm.res<-lm(ctmax ~ hi +STH.1 + STH.2 + STH.3 + STH.4 + STH.7 + STH.8 + STL.1 + 
                  STL.2 + STL.4 + STL.5 + STL.6 + STL.7 + STM.2 + STM.4 + STM.5 + 
                  STM.6 + STM.7,data=dp_noNA)
summary(null.lm.res)

## subtract the adjusted R^2 of the null model from the full model to get effect estimate
effect = 0.6417 - 0.6203 # 0.214 (2.14%)

### both QTL LM effect size estimate
## organize variables of interest into new data frame, remove missing data
dp_noNA<-data.frame(ctmax,hi,chr22_geno,chr15_geno,STH.1,STH.2,STH.3,STH.4,STH.7,STH.8,STL.1,STL.2,STL.4,STL.5,STL.6,STL.7,STM.2,STM.4,STM.5,STM.6,STM.7)
dp_noNA[dp_noNA == '-']<-NA
dp_noNA<-na.omit(dp_noNA)
dim(dp_noNA)
# full model adjusted R^2: 0.7126 (pval < 2e-16)
lm.res<-lm(ctmax ~ chr22_geno + chr15_geno + chr22_geno:chr15_geno + hi + STH.1 + STH.2 + STH.3 + STH.4 + STH.7 + STH.8 + STL.1 + 
             STL.2 + STL.4 + STL.5 + STL.6 + STL.7 + STM.2 + STM.4 + STM.5 + 
             STM.6 + STM.7,data=dp_noNA)
summary(lm.res)

## get the adjusted R^2 from the null model (i.e. without genotype as a coefficient)
# null adj R^2: 0.6306 (pval < 2e-16)
null.lm.res<-lm(ctmax ~ hi + STH.1 + STH.2 + STH.3 + STH.4 + STH.7 + STH.8 + STL.1 + 
                  STL.2 + STL.4 + STL.5 + STL.6 + STL.7 + STM.2 + STM.4 + STM.5 + 
                  STM.6 + STM.7,data=dp_noNA)
summary(null.lm.res)

## subtract the adjusted R^2 of the null model from the full model to get effect estimate
effect = 0.7126 - 0.6306 # 0.082 (8.2%)


### chr22:15 QTL interaction LM effect size estimate
## organize variables of interest into new data frame, remove missing data
library(mltools)
library(data.table)
library(rstatix)
dp_noNA<-data.frame(ctmax,hi,chr22_geno,chr15_geno,STH.1,STH.2,STH.3,STH.4,STH.7,STH.8,STL.1,STL.2,STL.4,STL.5,STL.6,STL.7,STM.2,STM.4,STM.5,STM.6,STM.7)
dp_noNA[dp_noNA == '-']<-NA
dp_noNA<-na.omit(dp_noNA)
dim(dp_noNA)

# full model adjusted R^2: 0.7126 (pval < 2e-16)
lm.res<-lm(ctmax ~ chr22_geno + chr15_geno + chr22_geno:chr15_geno + hi + STH.1 + STH.2 + STH.3 + STH.4 + STH.7 + STH.8 + STL.1 + 
             STL.2 + STL.4 + STL.5 + STL.6 + STL.7 + STM.2 + STM.4 + STM.5 + 
             STM.6 + STM.7,data=dp_noNA)
summary(lm.res)

## get the adjusted R^2 from the null model (i.e. without genotype as a coefficient)
# null adj R^2: 0.687 (pval < 2e-16)
null.lm.res<-lm(ctmax ~ chr22_geno + chr15_geno + hi + STH.1 + STH.2 + STH.3 + STH.4 + STH.7 + STH.8 + STL.1 + 
                  STL.2 + STL.4 + STL.5 + STL.6 + STL.7 + STM.2 + STM.4 + STM.5 + 
                  STM.6 + STM.7,data=dp_noNA)
summary(null.lm.res)

## subtract the adjusted R^2 of the null model from the full model to get effect estimate
effect = 0.7126 - 0.687 # 0.0256 (2.6%)

## eta squared to get effect size per variable
eta_squared(aov(ctmax ~ chr22_geno + chr15_geno + chr22_geno:chr15_geno + hi + STH.1 + STH.2 + STH.3 + STH.4 + STH.7 + STH.8 + STL.1 + STL.2 + STL.4 + STL.5 + STL.6 + STL.7 + STM.2 + STM.4 + STM.5 + STM.6 + STM.7, dp_noNA))


### chr22:chr15 QTL comparison of means w anova/tukey
lm.res<-lm(ctmax ~ chr22_geno + chr15_geno + chr22_geno:chr15_geno + hi + STH.1 + STH.2 + STH.3 + STH.4 + STH.7 + STH.8 + STL.1 + 
             STL.2 + STL.4 + STL.5 + STL.6 + STL.7 + STM.2 + STM.4 + STM.5 + 
             STM.6 + STM.7,data=dp_noNA)
anov <- aov(lm.res)
TukeyHSD(x=anov, 'chr22_geno', conf.level=0.95)
TukeyHSD(x=anov, 'chr15_geno', conf.level=0.95)
TukeyHSD(x=anov, 'chr22_geno:chr15_geno', conf.level=0.95)


###### CREATE A NICE MANHATTAN PLOT OF GENOME-WIDE LODs ######

### To generate an organized genome-wide manhattan plot, convert chromosome names to numbers using 
#     ./replace-chrom-name.py ./match_xmac_chroms_xbir-10x_chroms.txt <infile>
#   in ./Scripts/
scanone.hk.df_rename <- read.csv("ltreb-only-CTmax_17-1hot-site-tanks.scanone-hk_clean.tsv_chr-renamed.csv")

## load the qqman manhattan plot function
source("~/Box/Schumer_lab_resources/Project_files/Thermal_tolerance_projects/Scripts/adapted_qqman.R")

## Function: prepares data frame with chr,pos,lod columns to create manhattan plot
#  make sure 'cutoff_lod' is set to the LOD threshold
mh_plot<-function(chr_pos_lod){
  cpl.df <- as.data.frame(chr_pos_lod)
  data_trim<-{}
  # if using non-integer chrom names (i.e. not group numbers),
  # use this instead of following line
  # data_trim$CHR<-chr_to_group(as.character.factor(cpl.df$chr))
  data_trim$CHR<-cpl.df$chr
  data_trim$BP<-(cpl.df$pos)*1e6
  data_trim$P<-cpl.df$lod
  data_trim<-as.data.frame(data_trim)
  data_trim<-na.omit(data_trim)
  # data_trim$CHR <- as.character.factor(data_trim$CHR)
  data_trim$CHR<-data_trim$CHR
  manhattan(data_trim,genomewideline=cutoff_lod,suggestiveline=-1, cex.axis=1.5, cex.lab=4)
  
  return(data_trim)
}

## look at all of the peaks that surpass the genome-wide LOD threshold
scanone.hk.df_rename[scanone.hk.df_rename$lod>cutoff_lod,]
# ScyDAA6-2113-HRSCAF-2539:8811359  22 8.811359 4.366771
# ScyDAA6-2113-HRSCAF-2539:8831590  22 8.831590 4.361586
# ScyDAA6-2113-HRSCAF-2539:8851788  22 8.851788 4.356347

#  and all those under the 2-LOD interval
scanone.hk.df_rename[scanone.hk.df_rename$lod>lint2 & scanone.hk.df_rename$chr == "22" ,]
write.table(scanone.hk.df_rename[scanone.hk.df_rename$lod>lint2 & scanone.hk.df_rename$chr == "22" ,],'ltreb-ctmax_grp22-all-markers-under-2LOD_chr-pos-lod.tsv',sep="\t",row.names=F,quote=F)

## plot
pdf("chr22-ctmax-qtl_first-scan_genome-wide-manhattan.pdf",10.5,5.5)
mh_plot(scanone.hk.df_rename)                                         # genome-wide
title(main="LTREB CTmax QTL: QTL map ctmax~geno+hi+17sig.site.tanks")
dev.off()
pdf("chr22-ctmax-qtl_first-scan_chr22-manhattan.pdf",5.3,5.3)
mh_plot(scanone.hk.df_rename[scanone.hk.df_rename$chr == "22" ,])     # particular chromosome
dev.off()


###### CALCULATE CORRELATION BETWEEN CTMAX AND GENOME-WIDE ANCESTRY ######

## load phenotype data (individuals with excess hybrid_index and heterozygosity have already been removed)
pheno_data <- read.table('LTREB-only-CTmax-phenos-w-hi-het_1hot-site-tank_clean.tsv',header=T)
pheno_data <- pheno_data[pheno_data$hybrid_index>0.15 & pheno_data$hybrid_index<0.85,] # only keep individuals with hybrid index between 0.15 and 0.85
hi<-pheno_data$hybrid_index     # genome-wide X.malinche ancestry proportion
het<-pheno_data$heterozygosity  # genome-wide heterozygosity proportion
ctmax<-pheno_data$ctmax         # CTmax

## look at correlation between CTmax and genome-wide X. malinche ancestry (hybrid index)
model <- lm(ctmax~hi)
summary(model)
cor.test(ctmax,hi,method="pearson") # ctmax,hi: r=+0.1185137, pval=0.1472
# plot correlation 
pdf("CTmax-vs-genome-wide-ancestry_corr-plot.pdf",6,6)
plot(hi,ctmax,pch = 16, cex = 1.3, col = "light grey", main = "CTmax and genome-wide X. malinche ancestry",xlab="hybrid index",ylab="CTmax (Celsius)")
abline(34.7295,2.1980, lwd = 3, col="dark blue")
dev.off()

## look at correlation between CTmax and genome-wide heterozygosity)
model <- lm(ctmax~het)
summary(model)
cor.test(ctmax,het,method="pearson") # ctmax,het: r=-0.02285036, pval=0.7806
# plot correlation
pdf("CTmax-vs-genome-wide-heterozygosity_corr-plot.pdf",6,6)
plot(het,ctmax,pch = 16, cex = 1.3, col = "light salmon", main = "CTmax and genome-wide heterozygosity",xlab="heterozygosity",ylab="CTmax (Celsius)")
abline(35.9760,-0.2317, lwd = 3, col="dark blue")
dev.off()
