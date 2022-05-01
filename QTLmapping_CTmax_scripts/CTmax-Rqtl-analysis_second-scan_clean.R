### CTmax-Rqtl-analysis_second-scan_clean.R
#
# Script used to map second CTmax QTL with R/qtl in Payne et al 2021
# See Figure 2
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
#
#   Has an additional phenotype: first-scan chr22 QTL peak as variable

## load library
library(qtl)

###### LOAD DATA ######

## load 2nd QTL scan rqtl infile (csv) = 30244 (20k bp thinned) markers, 152 individuals, 31 phenotypes
data<-read.cross("csv",dir="Scripts/input_files","LTREB-CTmax-oct2020-ONLY_20k-thinned-genos_no-xtrm-hi_22onehot-site-tanks_ScyDAA6-2113-HRSCAF-2539.8811359-onehot.rqtl.csv", na.strings=c("NA"),estimate.map=FALSE)

# check that cross type, # individuals, and # markers are correct
# 30244 (20k bp thinned) markers, 152 individuals, 35 phenotypes
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
data_sub.ds <- readRDS(file = "./ltreb-qtl-data_sub.ds-est.rf_second-scan-w-qtl-covariate.rds")

## calculate conditional genotype probabilities given multipoint marker data
data_prob <- calc.genoprob(data_sub.ds)

write.cross(data_prob,file="ltreb-ctmax-qtl-filtered_22-1hot-site-tank_data-prob_second-scan_clean.csv")

saveRDS(data_prob,file="ltreb-ctmax-qtl-filtered_22-1hot-site-tank_data-prob_second-scan_clean.rds")


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

# 2nd scan
#geno <- as.factor(pull.geno(data_prob, "ScyDAA6-2113-HRSCAF-2539:8811359"))
#summary(lm(ctmax ~ geno + hi + STH.1 + STH.2 + STH.3 + STH.4 + STH.7 + STH.8 + STL.1 +
#             STL.2 + STL.4 + STL.5 + STL.6 + STL.7 + STM.2 + STM.4 + STM.5 +
#             STM.6 + STM.7))
# pull selected covariates (plus hybrid_index), including the interacting covariates
#  ctmax ~ hi + STH.1 + STH.2 + STH.3 + STH.4 + STH.7 + STH.8 + STL.1 +
#  STL.2 + STL.4 + STL.5 + STL.6 + STL.7 + STM.2 + STM.4 + STM.5 +
#  STM.6 + STM.7
covariates <- pull.pheno(data_prob, c("hybrid_index","STH.1","STH.2","STH.3","STH.4","STH.7","STH.8","STL.1","STL.2","STL.4","STL.5","STL.6","STL.7","STM.2",
                                      "STM.4","STM.5","STM.6","STM.7","homo_birch0.0","het1.0","homo_mal2.0"))
# add one-hot encoded genotype states as interaction terms
interactor <- pull.pheno(data_prob, c("homo_birch0.0","het1.0","homo_mal2.0"))

###### RUN SINGLE-QTL SCAN ######

## perform a single-qtl scan, with chr22 qtl genotypes as interaction covariate
scanone.hk <- scanone(data_prob, method="hk",addcovar=covariates,intcovar=interactor)
plot(scanone.hk,col="black")
summary(scanone.hk)
write.table(scanone.hk,file="ltreb-only-CTmax_17-1hot-site-tanks.scanone-hk_second-scan_clean.tsv",sep="\t",quote=FALSE)
#write.table(scanone.hk,file="ltreb-only-CTmax_22-1hot-site-tanks.scanone-hk_second-scan-w-qtl-covariate.tsv",sep="\t",quote=FALSE)
saveRDS(scanone.hk, file="ltreb-only-CTmax_17-1hot-site-tanks.scanone-hk_second-scan_clean.rds")

## run permutations to define the LOD threshold of significance
perm.hk <- scanone(data_prob, method="hk", addcovar=covariates, intcovar=interactor, n.perm=1000)
saveRDS(perm.hk, "LTREB-only-ctmax_perm_17-site-tank.hk_perm.hk_second-scan_clean.rds")

#perm.hk <- readRDS("LTREB-only-ctmax_perm_17-site-tank.hk_perm.hk_second-scan_clean.rds")
#quantile(perm.hk,c(0.8,0.9,0.95))

## get genome-wide LOD thresholds using genome-scan-adjusted p-values
summary(perm.hk,alpha=c(0.05, 0.1, 0.2))
## ctmax ~ chr22_qtl+hi+17sig.site.tanks
#lod
#5%  9.63
#10% 8.96
#20% 8.32
#1.5interval=8.96-1.5=7.46
cutoff_lod = 8.96

## look at highest LOD peak and the corresponding p-value on each chromosome
summary(scanone.hk, perms=perm.hk, pvalues=TRUE)

## look at peaks that surpass the 10% LOD threshold
summary(scanone.hk, perms=perm.hk, alpha=0.1, pvalues=TRUE)
summary(scanone.hk[scanone.hk$lod > cutoff_lod,])

# 2nd scan: chr15 ScyDAA6-5984-HRSCAF-6694:4408039 ScyDAA6-5984-HRSCAF-6694 4.41 9.1 0.088
# output all markers with LOD score above 10% LOD threshold
write.table(scanone.hk[scanone.hk$lod > cutoff_lod,],file="ltreb-only-CTmax_17-1hot-site-tanks.scanone-hk_peaks-above-0.1LODcutoff_second-scan_clean.tsv",sep="\t",quote=FALSE)

## plot chr22 peak
pdf("ltreb-only-CTmax_17-1hot-site-tanks_chr15-qtl_second-scan.pdf",8.4,8.9)
plot(scanone.hk,chr="ScyDAA6-5984-HRSCAF-6694") # chr15
abline(h=cutoff_lod, col="red", lwd=3)
title("15@4.41")
dev.off()


###### LOD peak support ######

## get the 1.5-LOD support interval for significant peaks
lint = cutoff_lod - 1.5   # 7.46
lodint(scanone.hk,"ScyDAA6-5984-HRSCAF-6694",drop=1.5,expandtomarkers = TRUE)
#ScyDAA6-5984-HRSCAF-6694:3880879 ScyDAA6-5984-HRSCAF-6694 3.880879 6.380639
#ScyDAA6-5984-HRSCAF-6694:4408039 ScyDAA6-5984-HRSCAF-6694 4.408039 9.101299
#ScyDAA6-5984-HRSCAF-6694:6054388 ScyDAA6-5984-HRSCAF-6694 6.054388 7.422079

## get the 2.0-LOD support interval for significant peaks
lint2 = cutoff_lod - 2   # 6.96
lodint(scanone.hk,"ScyDAA6-5984-HRSCAF-6694",drop=2,expandtomarkers = TRUE)
#ScyDAA6-5984-HRSCAF-6694:2958764 ScyDAA6-5984-HRSCAF-6694 2.958764 5.367273
#ScyDAA6-5984-HRSCAF-6694:4408039 ScyDAA6-5984-HRSCAF-6694 4.408039 9.101299
#ScyDAA6-5984-HRSCAF-6694:6156396 ScyDAA6-5984-HRSCAF-6694 6.156396 7.096667

## get 95% Bayes credible interval
bayesint(scanone.hk,"ScyDAA6-5984-HRSCAF-6694",0.95,expandtomarkers = TRUE)
#ScyDAA6-5984-HRSCAF-6694:2834112 ScyDAA6-5984-HRSCAF-6694 2.834112 6.373343
#ScyDAA6-5984-HRSCAF-6694:4408039 ScyDAA6-5984-HRSCAF-6694 4.408039 9.101299
#ScyDAA6-5984-HRSCAF-6694:6196895 ScyDAA6-5984-HRSCAF-6694 6.196895 6.657164


###### MAKE EFFECT PLOTS ######
## Load data_prob object
data_prob <- readRDS("ltreb-ctmax-qtl-filtered_22-1hot-site-tank_data-prob_second-scan_clean.rds")

## impute missing genotypes
data_sim <- sim.geno(data_prob,step=1,n.draws=16)

## plot estimated effect of QTL: phenotype vs marker genotype
malcol=rgb(0/255,0/255,175/255)
hetcol=rgb(100/255,0/255,175/255)
bircol=rgb(150/255,0/255,0/255)
effect_colors = c(malcol,hetcol,bircol)

# 1st scan
mar_chr22 <- find.marker(data_prob, chr='ScyDAA6-2113-HRSCAF-2539', pos=8.811359)
# simple effect plot
pdf("ltreb-only-CTmax_17-1hot-site-tanks_chr22-qtl-effect-plot.pdf",5.5,7)
effectplot(data_sim,mname1=mar_chr22,main='22@8.81',ylab = expression('CT'['max']*' ('*~degree*C*')'))
dev.off()
# effect plots with individuals as points
pdf("ltreb-only-CTmax_17-1hot-site-tanks_chr22-qtl-effect-plot-w-points.pdf",5.5,7)

plotPXG(data_prob,marker=mar_chr22,main='22@8.81',infer=FALSE, ylab = expression('CT'['max']*' ('*~degree*C*')'), col = effect_colors)

dev.off()

# 2nd scan
mar_chr15 <- find.marker(data_prob, chr='ScyDAA6-5984-HRSCAF-6694', pos=4.408039)
plotPXG(data_prob,marker=mar_chr15,main='15@4.41')
effectplot(data_sim,mname1=mar_chr15,main='15@4.41')

## plot the joint effects of the chr22 QTL and chr15 peak
pdf("ltreb-only-CTmax_17-1hot-site-tanks_chr22-qtl-chr15-peak-interaction-plot.pdf",5.5,7)
par(mfrow=c(1,1), las=1)
#effectplot(data_sim, mname1="ScyDAA6-2113-HRSCAF-2539@8.81", mname2="ScyDAA6-5984-HRSCAF-6694@4.41", main = "22@8.81 x 15@5.85")
effectplot(data_sim, mname2="ScyDAA6-2113-HRSCAF-2539@8.81", mname1="ScyDAA6-5984-HRSCAF-6694@4.41", geno1=c("MM","MB","BB"), geno2=c("MM","MB","BB"), var.flag = "group", main = "15@4.41 x 22@8.81", xlab = "Genotype at chr22 QTL", ylab = expression('CT'['max']*' ('*~degree*C*')'), col = effect_colors, add.legend = TRUE, legend.lab = "chr15 QTL")
dev.off()

pdf("ltreb-only-CTmax_17-1hot-site-tanks_chr22-qtl-chr15-peak-interaction-plot-w-points.pdf",9.8,8.2)
plotPXG(data_prob, marker=c(mar_chr22, mar_chr15))
dev.off()

## plot estimated effect of QTL: phenotype vs marker genotype
# 2nd scan
mar_chr15 <- find.marker(data_prob, chr='ScyDAA6-5984-HRSCAF-6694', pos=4.408039)
pdf("ltreb-only-CTmax_17-1hot-site-tanks_chr15-qtl-effect-plot_second-scan.pdf",8.4,8.9)
effectplot(data_sim,mname1=mar_chr15,main='15@4.41')
dev.off()
pdf("ltreb-only-CTmax_17-1hot-site-tanks_chr15-qtl-effect-plot-w-points_second-scan.pdf",8.4,8.9)
plotPXG(data_prob,marker=mar_chr15,main='15@4.41')
dev.off()

## get QTL effect on CTmax per genotype at peak
eff <- effectplot(data_sim,mname1=mar_chr15)
eff
# ScyDAA6-5984-HRSCAF-6694:4408039
#      mean      SE
#  AA: 35.81591, 0.2137864
#  AB: 35.64323, 0.1078122
#  BB: 36.27023, 0.1734309

## plot estimated QTL effects along chr15 (shows additive and dominance effects)
# effect summary for 15@4.41:
# a           d             se.a        se.d
# 0.22691415 -0.3949753     0.1383072 0.1772078
effect_table<-effectscan(data_sim, pheno.col="ctmax", c("ScyDAA6-5984-HRSCAF-6694"), draw=F, get.se=T)
write.table(effect_table,file="ltreb-only-CTmax_17-1hot-site-tanks.chr15-qtl_second-scan_effect-table.tsv",sep="\t",quote=FALSE)
# plot effect across chr22
pdf("ltreb-only-CTmax_17-1hot-site-tanks_chr22-additive-dominance-effects-plot.pdf",9.2,6.8)
effectscan(data_sim, pheno.col="ctmax", c("ScyDAA6-5984-HRSCAF-6694"), draw=T, get.se=T)
dev.off()


###### ROUGHLY ESTIMATE EFFECT OF QTL WITH AN LM ######

### get a very rough estimate of the QTL effect size with R^2
## load infile as data.frame and do correlations/lm with that
setwd("~/Box/Schumer_lab_resources/Project_files/Thermal_tolerance_projects/Data/LTREB_qtl/CTmax-QTL_second-scan_data")
data_prob.df <- read.csv("ltreb-ctmax-qtl-filtered_22-1hot-site-tank_data-prob_second-scan_clean.csv")

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
chr22_interactor <- data_prob.df$ScyDAA6.2113.HRSCAF.2539.8811359

# grab the genotypes of the locus at the QTL peak
chr15_geno <- data_prob.df$ScyDAA6.5984.HRSCAF.6694.4408039

## output genotypes for markers that fall under the QTL / within the 2-LOD interval: 2958764 - 6156396
firstcol = which(colnames(data_prob.df)=="ScyDAA6.5984.HRSCAF.6694.2958764")
lastcol = which(colnames(data_prob.df)=="ScyDAA6.5984.HRSCAF.6694.6156396")
qtl_genos <- data_prob.df[c(firstcol:lastcol)]
qtl_geno_outdf <- data.frame(ctmax,hi,STH.1,STH.2,STH.3,STH.4,STH.7,STH.8,STL.1,STL.2,STL.4,STL.5,STL.6,STL.7,STM.2,STM.4,STM.5,STM.6,STM.7,chr22_interactor,qtl_genos)
write.csv(qtl_geno_outdf,'./ltreb-ctmax-qtl_genos-all-markers-under-2LOD-chr15-qtl_second-scan.csv',row.names=F,quote=F)

## organize variables of interest into new data frame, remove missing data
dp_noNA<-data.frame(ctmax,hi,chr15_geno,STH.1,STH.2,STH.3,STH.4,STH.7,STH.8,STL.1,STL.2,STL.4,STL.5,STL.6,STL.7,STM.2,STM.4,STM.5,STM.6,STM.7,chr22_interactor)
dp_noNA[dp_noNA == '-']<-NA
dp_noNA<-na.omit(dp_noNA)
dim(dp_noNA)

## perform a linear regression with the full model and get the R^2 value
# full model adjusted R^2: 0.7126 (pval < 2e-16)
lm.res<-lm(ctmax ~ as.factor(chr15_geno) + hi + STH.1 + STH.2 + STH.3 + STH.4 + STH.7 + STH.8 + STL.1 +
             STL.2 + STL.4 + STL.5 + STL.6 + STL.7 + STM.2 + STM.4 + STM.5 +
             STM.6 + STM.7 + as.factor(chr22_interactor) + as.factor(chr15_geno):as.factor(chr22_interactor),data=dp_noNA)
summary(lm.res)
# an ANOVA should give the same result
aov.res <- aov(ctmax ~ chr15_geno + hi + STH.1 + STH.2 + STH.3 + STH.4 + STH.7 + STH.8 + STL.1 +
                 STL.2 + STL.4 + STL.5 + STL.6 + STL.7 + STM.2 + STM.4 + STM.5 +
                 STM.6 + STM.7 + chr22_interactor + chr15_geno:chr22_interactor,data=dp_noNA)
summary(aov.res)

## get the adjusted R^2 from the null model (i.e. without genotype as a coefficient)
# null adj R^2: 0.6666 (pval < 2e-16)
null.lm.res<-lm(ctmax ~ hi + STH.1 + STH.2 + STH.3 + STH.4 + STH.7 + STH.8 + STL.1 +
                  STL.2 + STL.4 + STL.5 + STL.6 + STL.7 + STM.2 + STM.4 + STM.5 +
                  STM.6 + STM.7 + as.factor(chr22_interactor),data=dp_noNA)
summary(null.lm.res)

## subtract the adjusted R^2 of the null model from the full model to get effect estimate
effect = 0.7126 - 0.6666 # 0.046 (4.6%)


## get log likelihood difference between focal and null models
summary(lm.res)
model2<-glm(as.numeric(ctmax)~as.numeric(hi)+
              as.factor(STH.1) + as.factor(STH.2) + as.factor(STH.3) +
              as.factor(STH.4) + as.factor(STH.7) + as.factor(STH.8) +
              as.factor(STL.1) + as.factor(STL.2) + as.factor(STL.4) +
              as.factor(STL.5) + as.factor(STL.6) + as.factor(STL.7) +
              as.factor(STM.2) + as.factor(STM.4) + as.factor(STM.5) +
              as.factor(STM.6) + as.factor(STM.7) + as.factor(chr22_interactor),data=dp_noNA,family="gaussian")
null<-logLik(model2)[1]

model1<-glm(as.numeric(ctmax)~as.numeric(hi)+
              as.factor(STH.1) + as.factor(STH.2) + as.factor(STH.3) +
              as.factor(STH.4) + as.factor(STH.7) + as.factor(STH.8) +
              as.factor(STL.1) + as.factor(STL.2) + as.factor(STL.4) +
              as.factor(STL.5) + as.factor(STL.6) + as.factor(STL.7) +
              as.factor(STM.2) + as.factor(STM.4) + as.factor(STM.5) +
              as.factor(STM.6) + as.factor(STM.7) + as.factor(chr22_interactor) +
              as.factor(chr22_interactor):as.factor(chr15_geno) + as.factor(chr15_geno),data=dp_noNA,family="gaussian")
focal<-logLik(model1)[1]

# calculate log likelihood difference
like_diff<-focal-null
like_diff # 13.5293


###### CREATE A NICE MANHATTAN PLOT OF GENOME-WIDE LODs ######

### To generate an organized genome-wide manhattan plot, convert chromosome names to numbers using
#     ./replace-chrom-name.py ./match_xmac_chroms_xbir-10x_chroms.txt <infile>
scanone.hk.df_rename <- read.csv("ltreb-only-CTmax_17-1hot-site-tanks.scanone-hk_second-scan_clean.tsv_chr-renamed.csv")

## load the qqman manhattan plot function
#source("~/Box/Schumer_lab_resources/Project_files/Thermal_tolerance_projects/Scripts/adapted_qqman.R")

library("qqman")

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
  data_trim$P<- (10^(-1*(cpl.df$lod)))
  data_trim$SNP<-cpl.df$snp
  data_trim<-as.data.frame(data_trim)
  data_trim<-na.omit(data_trim)
  # data_trim$CHR <- as.character.factor(data_trim$CHR)
  manhattan(data_trim, suggestiveline=NULL, genomewideline=NULL, ylim=c(0,10), cex.axis=1.5, cex.lab=1.5)

  return(data_trim)
}

mh_plot(scanone.hk.df_rename)

# significance thresholds
#5%  9.63
#10% 8.96
#20% 8.32
cutoff_lod = 8.96
abline(h=cutoff_lod, col="red")
abline(h=9.63, col="grey")

## look at all of the peaks that surpass the genome-wide LOD threshold
scanone.hk.df_rename[scanone.hk.df_rename$lod>cutoff_lod,]
#15 4.408039 9.101299 15:4.408039

lint1.5 = cutoff_lod - 1.5
#  and all those under the 1.5-LOD interval
scanone.hk.df_rename[scanone.hk.df_rename$lod>lint1.5 & scanone.hk.df_rename$chr == "15" ,]
write.table(scanone.hk.df_rename[scanone.hk.df_rename$lod>lint1.5 & scanone.hk.df_rename$chr == "15" ,],'ltreb-ctmax_grp15-all-markers-under-1.5LOD_chr-pos-lod.tsv',sep="\t",row.names=F,quote=F)

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
