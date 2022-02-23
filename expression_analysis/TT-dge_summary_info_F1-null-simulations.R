## Summary of brain and liver differential expression results (from DESeq2 anlaysis)
## for simulated null F1 expression data
## adapted from ./TT-dge_summary_info.R
## FDR-adjusted p-value cutoff = 0.1
## Payne et al 2021

### BRAIN RNAseq (X. birchmanni pseudoalignment)

brain_dge_og <- read.csv("Scripts/input_files/nullF1_TT-brain-xbirch-gtf_DGE_lfc-shr_all.csv")
# first subset by genes that are expressed in both bir and mal under both thermal contexts
brain_dge <- subset(brain_dge_og,!is.na(padj_res.22CbirVmal) & !is.na(padj_res.33CbirVmal))
dim(brain_dge) #22931

## BRAIN 22C

## brain 22c: how many genes that are sig diff expressed between xbirVxmal have xbir-like f1 expression?
brain_dge_22c_f1bir <- subset(brain_dge, padj_res.22CbirVmal < 0.1 & !(padj_res.22CbirVF1 < 0.1) & padj_res.22CF1Vmal < 0.1 ) # f1 ge bir-like = 400
dim(subset(brain_dge_22c_f1bir,LFC_res.22CbirVmal < 0 & LFC_res.22CF1Vmal < 0)) # 52 xbir-like f1s have sig lower exp compared to xmal  [f1&bir < mal]
dim(subset(brain_dge_22c_f1bir,LFC_res.22CbirVmal > 0 & LFC_res.22CF1Vmal > 0)) # 348 xbir-like f1s have sig higher exp compared to xmal [f1&bir > mal]
# total: 52+348=400

## brain 22c: how many genes that are sig diff expressed between xbirVxmal have xmal-like f1 expression?
brain_dge_22c_f1mal <- subset(brain_dge, padj_res.22CbirVmal < 0.1 & padj_res.22CbirVF1 < 0.1 & !(padj_res.22CF1Vmal < 0.1) ) # f1 ge mal-like = 310
dim(subset(brain_dge_22c_f1mal,LFC_res.22CbirVmal < 0 & LFC_res.22CbirVF1 < 0)) # 170 xmal-like f1s have sig higher exp compared to xbir [f1&mal > bir]
dim(subset(brain_dge_22c_f1mal,LFC_res.22CbirVmal > 0 & LFC_res.22CbirVF1 > 0)) # 140 xmal-like f1s have sig lower exp compared to xbir  [f1&mal < bir]
# total: 140+170=310

## brain 22c: how many genes that are sig diff expressed for all comparisons have:
brain_dge_22c_allsig <- subset(brain_dge, padj_res.22CbirVmal < 0.1 & padj_res.22CbirVF1 < 0.1 & padj_res.22CF1Vmal < 0.1 ) # 22c all groups sig = 58
#intermediate f1 expression? 57
dim(subset(brain_dge_22c_allsig,LFC_res.22CbirVmal > 0 & LFC_res.22CbirVF1 > 0 & LFC_res.22CF1Vmal > 0 )) # bir > f1 > mal = 23
dim(subset(brain_dge_22c_allsig,LFC_res.22CbirVmal < 0 & LFC_res.22CbirVF1 < 0 & LFC_res.22CF1Vmal < 0 )) # bir < f1 < mal = 34
#sig high f1 misexpression? 0
dim(subset(brain_dge_22c_allsig, LFC_res.22CbirVF1 < 0 & LFC_res.22CF1Vmal > 0 )) # bir&mal < f1 = 0
#sig low f1 misexpression? 1
dim(subset(brain_dge_22c_allsig, LFC_res.22CbirVF1 > 0 & LFC_res.22CF1Vmal < 0 )) # bir&mal > f1 = 1

## brain 22c: how many genes are misexpressed in f1s? 231
f1_22c_sig_brain <- subset(brain_dge, padj_res.22CbirVF1 < 0.1 & padj_res.22CF1Vmal < 0.1 ) # 231
dim(subset(f1_22c_sig_brain, LFC_res.22CbirVF1 < 0 & LFC_res.22CF1Vmal > 0 )) # bir&mal < f1 = 72
dim(subset(f1_22c_sig_brain, LFC_res.22CbirVF1 > 0 & LFC_res.22CF1Vmal < 0 )) # bir&mal > f1 = 159


### BRAIN 33C

## brain 33c: how many genes that are sig diff expressed between xbirVxmal have xbir-like f1 expression?
brain_dge_33c_f1bir <- subset(brain_dge, padj_res.33CbirVmal < 0.1 & !(padj_res.33CbirVF1 < 0.1) & padj_res.33CF1Vmal < 0.1 ) # f1 ge bir-like = 180
dim(subset(brain_dge_33c_f1bir,LFC_res.33CbirVmal < 0 & LFC_res.33CF1Vmal < 0)) # 31 xbir-like f1s have sig lower exp compared to xmal  [f1&bir < mal]
dim(subset(brain_dge_33c_f1bir,LFC_res.33CbirVmal > 0 & LFC_res.33CF1Vmal > 0)) # 149 xbir-like f1s have sig higher exp compared to xmal [f1&bir > mal]
# total: 31+149=180

## brain 33c: how many genes that are sig diff expressed between xbirVxmal have xmal-like f1 expression?
brain_dge_33c_f1mal <- subset(brain_dge, padj_res.33CbirVmal < 0.1 & padj_res.33CbirVF1 < 0.1 & !(padj_res.33CF1Vmal < 0.1) ) # f1 ge mal-like = 294
dim(subset(brain_dge_33c_f1mal,LFC_res.33CbirVmal < 0 & LFC_res.33CbirVF1 < 0)) # 237 xmal-like f1s have sig higher exp compared to xbir [f1&mal > bir]
dim(subset(brain_dge_33c_f1mal,LFC_res.33CbirVmal > 0 & LFC_res.33CbirVF1 > 0)) # 57 xmal-like f1s have sig lower exp compared to xbir  [f1&mal < bir]
# total: 237+57=294

## brain 33c: how many genes that are sig diff expressed for all comparisons have:
brain_dge_33c_allsig <- subset(brain_dge, padj_res.33CbirVmal < 0.1 & padj_res.33CbirVF1 < 0.1 & padj_res.33CF1Vmal < 0.1 ) # 33c all groups sig = 50
#intermediate f1 expression? 49
dim(subset(brain_dge_33c_allsig,LFC_res.33CbirVmal > 0 & LFC_res.33CbirVF1 > 0 & LFC_res.33CF1Vmal > 0 )) # bir > f1 > mal = 35
dim(subset(brain_dge_33c_allsig,LFC_res.33CbirVmal < 0 & LFC_res.33CbirVF1 < 0 & LFC_res.33CF1Vmal < 0 )) # bir < f1 < mal = 14
#sig high f1 misexpression? 0
dim(subset(brain_dge_33c_allsig, LFC_res.33CbirVF1 < 0 & LFC_res.33CF1Vmal > 0 )) # bir&mal < f1 = 0
#sig low f1 misexpression? 1
dim(subset(brain_dge_33c_allsig, LFC_res.33CbirVF1 > 0 & LFC_res.33CF1Vmal < 0 )) # bir&mal > f1 = 1

## brain 33c: how many genes are misexpressed in f1s? 79
f1_33c_sig_brain <- subset(brain_dge, padj_res.33CbirVF1 < 0.1 & padj_res.33CF1Vmal < 0.1 )
dim(subset(f1_33c_sig_brain, LFC_res.33CbirVF1 < 0 & LFC_res.33CF1Vmal > 0 )) # bir&mal < f1 = 14
dim(subset(f1_33c_sig_brain, LFC_res.33CbirVF1 > 0 & LFC_res.33CF1Vmal < 0 )) # bir&mal > f1 = 65

## null brain misexpression summary
# total genes: 24173
# 22C: 231/22931 = 0.0101 = 1.01%
# 33C: 79/22931 = 0.00345 = 0.345%
# 306/22931= 1.33% genes are misexpressed in f1s under at least one thermal context
# 5/355 = 1.4% genes that were responsive to temperature in mal and/or birch were misexpressed in F1s under at least one thermal context

## how many of the genes with 33cV22c dge in Xbirch and/or Xmal are misexpressed in F1s under at least one thermal context
# intersect of genes dge between 22v33c birch and 22v33c mal
# brain
sig22v33_bir_brain <- subset(brain_dge, padj_res.bir33cV22c < 0.1) #285
sig22v33_mal_brain <- subset(brain_dge, padj_res.mal33cV22c < 0.1) #154
intersect(sig22v33_bir_brain$Gene,sig22v33_mal_brain$Gene) #84
dge33v22_xbir_xmal_brain <- union(sig22v33_bir_brain$Gene,sig22v33_mal_brain$Gene) #355
#285+154-84=355 genes are dge in birch and/or mal between 22c vs 33c

# pull all misexpressed genes from both thermal contexts
f1_22c_brain_misexp <- subset(f1_22c_sig_brain, (LFC_res.22CbirVF1 < 0 & LFC_res.22CF1Vmal > 0) | (LFC_res.22CbirVF1 > 0 & LFC_res.22CF1Vmal < 0 )) #231
f1_33c_brain_misexp <- subset(f1_33c_sig_brain, (LFC_res.33CbirVF1 < 0 & LFC_res.33CF1Vmal > 0) | (LFC_res.33CbirVF1 > 0 & LFC_res.33CF1Vmal < 0 )) #79
intersect(f1_22c_brain_misexp$Gene,f1_33c_brain_misexp$Gene) #4
f1misexp_brain <- union(f1_22c_brain_misexp$Gene,f1_33c_brain_misexp$Gene) #306
#231+79-4=306 (306/22931= 1.33%) genes are misexpressed in f1s under either 22c, 33c, or both

# get the intersect of all 33v22c genes in xbir/xmal (355) and the misexpressed set
intersect(dge33v22_xbir_xmal_brain,f1misexp_brain) #5

# overall, 5/355 (1.4%) genes that were responsive to temperature in xbirch and/or xmal were misexpressed in f1s under at least one thermal context



### LIVER RNAseq (X. birchmanni pseudoalignment)

liver_dge_og <- read.csv("Scripts/input_files/nullF1_TT-liver-xbirch-gtf_DGE_lfc-shr_all.csv")
liver_dge <- subset(liver_dge_og,!is.na(padj_res.22CbirVmal) & !is.na(padj_res.33CbirVmal))
dim(liver_dge) #17525

### LIVER 22c

## liver 22c: how many genes that are sig diff expressed between xbirVxmal have xbir-like f1 expression?
liver_dge_22c_f1bir <- subset(liver_dge, padj_res.22CbirVmal < 0.1 & !(padj_res.22CbirVF1 < 0.1) & padj_res.22CF1Vmal < 0.1 ) # f1 ge bir-like = 318
dim(subset(liver_dge_22c_f1bir,LFC_res.22CbirVmal < 0 & LFC_res.22CF1Vmal < 0)) # 11 xbir-like f1s have sig lower exp compared to xmal  [f1&bir < mal]
dim(subset(liver_dge_22c_f1bir,LFC_res.22CbirVmal > 0 & LFC_res.22CF1Vmal > 0)) # 307 xbir-like f1s have sig higher exp compared to xmal [f1&bir > mal]
# total: 11+307=318

## liver 22c: how many genes that are sig diff expressed between xbirVxmal have xbir-like f1 expression?
liver_dge_22c_f1mal <- subset(liver_dge, padj_res.22CbirVmal < 0.1 & padj_res.22CbirVF1 < 0.1 & !(padj_res.22CF1Vmal < 0.1) ) # f1 ge mal-like = 378
dim(subset(liver_dge_22c_f1mal,LFC_res.22CbirVmal < 0 & LFC_res.22CbirVF1 < 0)) # 292 xmal-like f1s have sig higher exp compared to xbir [f1&mal > bir]
dim(subset(liver_dge_22c_f1mal,LFC_res.22CbirVmal > 0 & LFC_res.22CbirVF1 > 0)) # 86 xmal-like f1s have sig lower exp compared to xbir  [f1&mal < bir]
# total: 292+86=378

## liver 22c: how many genes that are sig diff expressed for all comparisons have:
liver_dge_22c_allsig <- subset(liver_dge, padj_res.22CbirVmal < 0.1 & padj_res.22CbirVF1 < 0.1 & padj_res.22CF1Vmal < 0.1 ) # 22c all groups sig = 34
#intermediate f1 expression? 34
dim(subset(liver_dge_22c_allsig,LFC_res.22CbirVmal > 0 & LFC_res.22CbirVF1 > 0 & LFC_res.22CF1Vmal > 0 )) # bir > f1 > mal = 25
dim(subset(liver_dge_22c_allsig,LFC_res.22CbirVmal < 0 & LFC_res.22CbirVF1 < 0 & LFC_res.22CF1Vmal < 0 )) # bir < f1 < mal = 9
#sig high f1 misexpression? 0
dim(subset(liver_dge_22c_allsig, LFC_res.22CbirVF1 < 0 & LFC_res.22CF1Vmal > 0 )) # bir&mal < f1 = 0
#sig low f1 misexpression? 0
dim(subset(liver_dge_22c_allsig, LFC_res.22CbirVF1 > 0 & LFC_res.22CF1Vmal < 0 )) # bir&mal > f1 = 0

## liver 22c: how many genes are misexpressed in f1s? 3 / 24173
f1_22c_sig_liver <- subset(liver_dge, padj_res.22CbirVF1 < 0.1 & padj_res.22CF1Vmal < 0.1 ) #37
dim(subset(f1_22c_sig_liver, LFC_res.22CbirVF1 < 0 & LFC_res.22CF1Vmal > 0 )) # bir&mal < f1 = 0
dim(subset(f1_22c_sig_liver, LFC_res.22CbirVF1 > 0 & LFC_res.22CF1Vmal < 0 )) # bir&mal > f1 = 3


### LIVER 33c

## liver 33c: how many genes that are sig diff expressed between xbirVxmal have xbir-like f1 expression?
liver_dge_33c_f1bir <- subset(liver_dge, padj_res.33CbirVmal < 0.1 & !(padj_res.33CbirVF1 < 0.1) & padj_res.33CF1Vmal < 0.1 ) # f1 ge bir-like = 65
dim(subset(liver_dge_33c_f1bir,LFC_res.33CbirVmal < 0 & LFC_res.33CF1Vmal < 0)) # 0 xbir-like f1s have sig lower exp compared to xmal  [f1&bir < mal]
dim(subset(liver_dge_33c_f1bir,LFC_res.33CbirVmal > 0 & LFC_res.33CF1Vmal > 0)) # 65 xbir-like f1s have sig higher exp compared to xmal [f1&bir > mal]
# total: 0+65=65

## liver 33c: how many genes that are sig diff expressed between xbirVxmal have xmal-like f1 expression?
liver_dge_33c_f1mal <- subset(liver_dge, padj_res.33CbirVmal < 0.1 & padj_res.33CbirVF1 < 0.1 & !(padj_res.33CF1Vmal < 0.1) ) # f1 ge mal-like = 829
dim(subset(liver_dge_33c_f1mal,LFC_res.33CbirVmal < 0 & LFC_res.33CbirVF1 < 0)) # 617 xmal-like f1s have sig higher exp compared to xbir [f1&mal > bir]
dim(subset(liver_dge_33c_f1mal,LFC_res.33CbirVmal > 0 & LFC_res.33CbirVF1 > 0)) # 212 xmal-like f1s have sig lower exp compared to xbir  [f1&mal < bir]
# total: 617+212=829

## liver 33c: how many genes that are sig diff expressed for all comparisons have:
liver_dge_33c_allsig <- subset(liver_dge, padj_res.33CbirVmal < 0.1 & padj_res.33CbirVF1 < 0.1 & padj_res.33CF1Vmal < 0.1 ) # 33c all groups sig = 38
#intermediate f1 expression? 38
dim(subset(liver_dge_33c_allsig,LFC_res.33CbirVmal > 0 & LFC_res.33CbirVF1 > 0 & LFC_res.33CF1Vmal > 0 )) # bir > f1 > mal = 38
dim(subset(liver_dge_33c_allsig,LFC_res.33CbirVmal < 0 & LFC_res.33CbirVF1 < 0 & LFC_res.33CF1Vmal < 0 )) # bir < f1 < mal = 0
#sig high f1 misexpression? 0
dim(subset(liver_dge_33c_allsig, LFC_res.33CbirVF1 < 0 & LFC_res.33CF1Vmal > 0 )) # bir&mal < f1 = 0
#sig low f1 misexpression? 0
dim(subset(liver_dge_33c_allsig, LFC_res.33CbirVF1 > 0 & LFC_res.33CF1Vmal < 0 )) # bir&mal > f1 = 0

## liver 33c: how many genes are misexpressed in f1s? 0 / 24173
f1_33c_sig_liver <- subset(liver_dge, padj_res.33CbirVF1 < 0.1 & padj_res.33CF1Vmal < 0.1 ) #38
dim(subset(f1_33c_sig_liver, LFC_res.33CbirVF1 < 0 & LFC_res.33CF1Vmal > 0 )) # bir&mal < f1 = 0
dim(subset(f1_33c_sig_liver, LFC_res.33CbirVF1 > 0 & LFC_res.33CF1Vmal < 0 )) # bir&mal > f1 = 0

## liver misexpression summary
# total genes: 24173
# 22C: 3/17525 = 0.00017 = 0.017%
# 33C: 0/17525 = 0.0000 = 0.00%
# 3/17525= 0.017% genes are misexpressed in f1s under at least one thermal context
# 0/699 = 0.00% genes that were responsive to temperature in mal and/or birch were misexpressed in F1s under at least one thermal context

## how many of the genes with 33cV22c dge in Xbirch and/or Xmal are misexpressed in F1s under at least one thermal context?
# intersect of genes dge between 22v33c brich and 22v33c mal
# liver

sig22v33_bir_liver <- subset(liver_dge, padj_res.bir33cV22c < 0.1) #644
sig22v33_mal_liver <- subset(liver_dge, padj_res.mal33cV22c < 0.1) #87
intersect(sig22v33_bir_liver$Gene,sig22v33_mal_liver$Gene) #32
dge33v22_xbir_xmal_liver <- union(sig22v33_bir_liver$Gene,sig22v33_mal_liver$Gene) #699
# 644+87-32=699 genes are dge in birch and/or mal between 22c vs 33c

# pull all misexpressed genes from both thermal contexts
f1_22c_liver_misexp <- subset(f1_22c_sig_liver, (LFC_res.22CbirVF1 < 0 & LFC_res.22CF1Vmal > 0) | (LFC_res.22CbirVF1 > 0 & LFC_res.22CF1Vmal < 0 )) #3
f1_33c_liver_misexp <- subset(f1_33c_sig_liver, (LFC_res.33CbirVF1 < 0 & LFC_res.33CF1Vmal > 0) | (LFC_res.33CbirVF1 > 0 & LFC_res.33CF1Vmal < 0 )) #0
intersect(f1_22c_liver_misexp$Gene,f1_33c_liver_misexp$Gene) #0
f1misexp_liver <- union(f1_22c_liver_misexp$Gene,f1_33c_liver_misexp$Gene) #3
#3+0-0=3 (0.01%) genes are misexpressed in f1s under either 22c, 33c, or both

# get the intersect of all 33v22c genes in xbir/xmal (134) and the misexxpressed set
intersect(dge33v22_xbir_xmal_liver,f1misexp_liver) #0

# overall, 0/699 (0.0%) genes that were responsive to temperature in xbirch and/or xmal were misexpressed in f1s under at least one thermal context

