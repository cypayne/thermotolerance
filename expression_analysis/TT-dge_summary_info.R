## TT-dge_summary_info.R
## 
## Summary of brain and liver differential expression results (from DESeq2 anlaysis)
## FDR-adjusted p-value cutoff = 0.1
##
## cyp I-2022 

### BRAIN RNAseq (X. birchmanni pseudoalignment)

### BRAIN 33v22c
## summary of all genes differentially expressed between temperature treatments
xbir_brain_dge_33v22c <- read.csv("input_files/TT-brain-33cV22c-xbirch-gtf_DGE_lfc-shr_all.csv")
xbir_brain_sig <- subset(xbir_brain_dge_33v22c,padj_res.bir33cV22c < 0.1)
f1_brain_sig <- subset(xbir_brain_dge_33v22c,padj_res.F133cV22c < 0.1)
xmal_brain_sig <- subset(xbir_brain_dge_33v22c,padj_res.mal33cV22c < 0.1)
dim(xbir_brain_sig) # 882
dim(f1_brain_sig)   # 2318
dim(xmal_brain_sig) # 979

# genes uniquely differentially expressed in response to thermal stress per group
all_sig <- subset(xbir_brain_dge_33v22c,padj_res.bir33cV22c < 0.1 & padj_res.F133cV22c < 0.1 & padj_res.mal33cV22c < 0.1)
bir_only_sig <- subset(xbir_brain_dge_33v22c,padj_res.bir33cV22c < 0.1 & !(padj_res.F133cV22c < 0.1) & !(padj_res.mal33cV22c < 0.1))
f1_only_sig <- subset(xbir_brain_dge_33v22c,!(padj_res.bir33cV22c < 0.1) & padj_res.F133cV22c < 0.1 & !(padj_res.mal33cV22c < 0.1))
mal_only_sig <- subset(xbir_brain_dge_33v22c,!(padj_res.bir33cV22c < 0.1) & !(padj_res.F133cV22c < 0.1) & padj_res.mal33cV22c < 0.1)
bir_mal_only_sig  <- subset(xbir_brain_dge_33v22c,padj_res.bir33cV22c < 0.1 & !(padj_res.F133cV22c < 0.1) & padj_res.mal33cV22c < 0.1)
bir_mal_sig  <- subset(xbir_brain_dge_33v22c,padj_res.bir33cV22c < 0.1 & padj_res.mal33cV22c < 0.1)
bir_mal_sig_negLFC <- subset(xbir_brain_dge_33v22c,padj_res.bir33cV22c < 0.1 & padj_res.mal33cV22c < 0.1 & LFC_res.bir33cV22c < 0 & LFC_res.mal33cV22c < 0)
bir_mal_sig_posLFC <- subset(xbir_brain_dge_33v22c,padj_res.bir33cV22c < 0.1 & padj_res.mal33cV22c < 0.1 & LFC_res.bir33cV22c > 0 & LFC_res.mal33cV22c > 0)
subset(xbir_brain_dge_33v22c,padj_res.bir33cV22c < 0.1 & padj_res.mal33cV22c < 0.1 & LFC_res.bir33cV22c < 0 & LFC_res.mal33cV22c > 0) #4
subset(xbir_brain_dge_33v22c,padj_res.bir33cV22c < 0.1 & padj_res.mal33cV22c < 0.1 & LFC_res.bir33cV22c > 0 & LFC_res.mal33cV22c < 0) #0

dim(all_sig) # 191
dim(bir_only_sig) # 329
dim(f1_only_sig)  # 1472
dim(mal_only_sig) # 479
dim(bir_mal_only_sig) # 167
dim(bir_mal_sig) # 358 = 354 of these are dge in the same direction, 4 are opposite. All 4 are downreg at 33c in bir but upreg at 33c in mal.
dim(bir_mal_sig_negLFC) # 94
dim(bir_mal_sig_posLFC) # 260


### BRAIN 22c
brain_dge_22c_og <- read.csv("Scripts/input_files/TT-brain-22c-xbirch-gtf_DGE_lfc-shr_all.csv_with-F1info.csv")
# first subset by genes that are expressed in both bir and mal under both thermal contexts
brain_dge_22c <- subset(brain_dge_22c_og,!is.na(padj_res.22CbirVmal))
dim(brain_dge_22c) # 22065

## brain 22c: how many genes that are sig diff expressed between xbirVxmal have xbir-like f1 expression?
brain_dge_22c_f1bir <- subset(brain_dge_22c, padj_res.22CbirVmal < 0.1 & !(padj_res.22CbirVF1 < 0.1) & padj_res.22CF1Vmal < 0.1 ) # f1 ge bir-like = 778
dim(subset(brain_dge_22c_f1bir,LFC_res.22CbirVmal < 0 & LFC_res.22CF1Vmal < 0)) # 280 xbir-like f1s have sig lower exp compared to xmal  [f1&bir < mal]
dim(subset(brain_dge_22c_f1bir,LFC_res.22CbirVmal > 0 & LFC_res.22CF1Vmal > 0)) # 498 xbir-like f1s have sig higher exp compared to xmal [f1&bir > mal]
# total: 280+498=778

## brain 22c: how many genes that are sig diff expressed between xbirVxmal have xbir-like f1 expression?
brain_dge_22c_f1mal <- subset(brain_dge_22c, padj_res.22CbirVmal < 0.1 & padj_res.22CbirVF1 < 0.1 & !(padj_res.22CF1Vmal < 0.1) ) # f1 ge mal-like = 843
dim(subset(brain_dge_22c_f1mal,LFC_res.22CbirVmal < 0 & LFC_res.22CbirVF1 < 0)) # 470 xmal-like f1s have sig higher exp compared to xbir [f1&mal > bir]
dim(subset(brain_dge_22c_f1mal,LFC_res.22CbirVmal > 0 & LFC_res.22CbirVF1 > 0)) # 373 xmal-like f1s have sig lower exp compared to xbir  [f1&mal < bir]
# total: 470+373=843

## brain 22c: how many genes that are sig diff expressed for all comparisons have:
brain_dge_22c_allsig <- subset(brain_dge_22c, padj_res.22CbirVmal < 0.1 & padj_res.22CbirVF1 < 0.1 & padj_res.22CF1Vmal < 0.1 ) # 22c all groups sig = 83
#intermediate f1 expression? 40
dim(subset(brain_dge_22c_allsig,LFC_res.22CbirVmal > 0 & LFC_res.22CbirVF1 > 0 & LFC_res.22CF1Vmal > 0 )) # bir > f1 > mal = 15
dim(subset(brain_dge_22c_allsig,LFC_res.22CbirVmal < 0 & LFC_res.22CbirVF1 < 0 & LFC_res.22CF1Vmal < 0 )) # bir < f1 < mal = 25
#sig high f1 misexpression? 33
dim(subset(brain_dge_22c_allsig, LFC_res.22CbirVF1 < 0 & LFC_res.22CF1Vmal > 0 )) # bir&mal < f1 = 33
#sig low f1 misexpression? 10
dim(subset(brain_dge_22c_allsig, LFC_res.22CbirVF1 > 0 & LFC_res.22CF1Vmal < 0 )) # bir&mal > f1 = 10

## brain 22c: how many genes are misexpressed in f1s? 521 / 22065
f1_sig <- subset(brain_dge_22c, padj_res.22CbirVF1 < 0.1 & padj_res.22CF1Vmal < 0.1 )
dim(subset(f1_sig, LFC_res.22CbirVF1 < 0 & LFC_res.22CF1Vmal > 0 )) # bir&mal < f1 = 246
dim(subset(f1_sig, LFC_res.22CbirVF1 > 0 & LFC_res.22CF1Vmal < 0 )) # bir&mal > f1 = 275


### BRAIN 33c
brain_dge_33c_og <- read.csv("Scripts/input_files/TT-brain-33c-xbirch-gtf_DGE_lfc-shr_all.csv_with-F1info.csv")
# first subset by genes that are expressed in both bir and mal under both thermal contexts
brain_dge_33c <- subset(brain_dge_33c_og,!is.na(padj_res.33CbirVmal))
dim(brain_dge_33c) # 22065

## brain 33c: how many genes that are sig diff expressed between xbirVxmal have xbir-like f1 expression?
brain_dge_33c_f1bir <- subset(brain_dge_33c, padj_res.33CbirVmal < 0.1 & !(padj_res.33CbirVF1 < 0.1) & padj_res.33CF1Vmal < 0.1 ) # f1 ge bir-like = 340
dim(subset(brain_dge_33c_f1bir,LFC_res.33CbirVmal < 0 & LFC_res.33CF1Vmal < 0)) # 137 xbir-like f1s have sig lower exp compared to xmal  [f1&bir < mal]
dim(subset(brain_dge_33c_f1bir,LFC_res.33CbirVmal > 0 & LFC_res.33CF1Vmal > 0)) # 203 xbir-like f1s have sig higher exp compared to xmal [f1&bir > mal]
# total: 137+203=340

## brain 33c: how many genes that are sig diff expressed between xbirVxmal have xbir-like f1 expression?
brain_dge_33c_f1mal <- subset(brain_dge_33c, padj_res.33CbirVmal < 0.1 & padj_res.33CbirVF1 < 0.1 & !(padj_res.33CF1Vmal < 0.1) ) # f1 ge mal-like = 1077
dim(subset(brain_dge_33c_f1mal,LFC_res.33CbirVmal < 0 & LFC_res.33CbirVF1 < 0)) # 567 xmal-like f1s have sig higher exp compared to xbir [f1&mal > bir]
dim(subset(brain_dge_33c_f1mal,LFC_res.33CbirVmal > 0 & LFC_res.33CbirVF1 > 0)) # 510 xmal-like f1s have sig lower exp compared to xbir  [f1&mal < bir]
# total: 567+510=1077

## brain 33c: how many genes that are sig diff expressed for all comparisons have:
brain_dge_33c_allsig <- subset(brain_dge_33c, padj_res.33CbirVmal < 0.1 & padj_res.33CbirVF1 < 0.1 & padj_res.33CF1Vmal < 0.1 ) # 33c all groups sig = 58
#intermediate f1 expression? 43
dim(subset(brain_dge_33c_allsig,LFC_res.33CbirVmal > 0 & LFC_res.33CbirVF1 > 0 & LFC_res.33CF1Vmal > 0 )) # bir > f1 > mal = 22
dim(subset(brain_dge_33c_allsig,LFC_res.33CbirVmal < 0 & LFC_res.33CbirVF1 < 0 & LFC_res.33CF1Vmal < 0 )) # bir < f1 < mal = 21
#sig high f1 misexpression? 33
dim(subset(brain_dge_33c_allsig, LFC_res.33CbirVF1 < 0 & LFC_res.33CF1Vmal > 0 )) # bir&mal < f1 = 10
#sig low f1 misexpression? 10
dim(subset(brain_dge_33c_allsig, LFC_res.33CbirVF1 > 0 & LFC_res.33CF1Vmal < 0 )) # bir&mal > f1 = 5

## brain 33c: how many genes are misexpressed in f1s? 187
f1_sig <- subset(brain_dge_33c, padj_res.33CbirVF1 < 0.1 & padj_res.33CF1Vmal < 0.1 )
dim(subset(f1_sig, LFC_res.33CbirVF1 < 0 & LFC_res.33CF1Vmal > 0 )) # bir&mal < f1 = 92
dim(subset(f1_sig, LFC_res.33CbirVF1 > 0 & LFC_res.33CF1Vmal < 0 )) # bir&mal > f1 = 52

## brain misexpression summary
# total genes: 22065
# 22C: 521/22065 = 0.0236 = 2.36%
# 33C: 144/22065 = 0.0065 = 0.65%
# 648/22065 = 2.94% genes are misexpressed in f1s under at least one thermal context
# 136/1503 = 9.0% genes that were responsive to temperature in mal and/or birch were misexpressed in F1s under at least one thermal context

## how many of the genes with 33cV22c dge in Xbirch and/or Xmal are misexpressed in F1s under at least one thermal context
# intersect of genes dge between 22v33c birch and 22v33c mal
# brain
# xbirch annot dge data
exp22v33c_brain <- read.csv("input_files/TT-brain-33cV22c-xbirch-gtf_DGE_lfc-shr_all.csv",header=T)

sig22v33_bir <- subset(exp22v33c_brain, padj_res.bir33cV22c < 0.1) #882
sig22v33_mal <- subset(exp22v33c_brain, padj_res.mal33cV22c < 0.1) #979
intersect(sig22v33_bir$Gene,sig22v33_mal$Gene) #358
dge33v22_xbir_xmal <- union(sig22v33_bir$Gene,sig22v33_mal$Gene)
#882+979-358=1503 genes are dge in birch and/or mal between 22c vs 33c

# pull all misexpressed genes from both thermal contexts
f1_22c_brain <- read.csv("input_files/TT-brain-22c-xbirch-gtf_DGE_lfc-shr_all.csv_with-F1info.csv",header=T)
f1_33c_brain <- read.csv("input_files/TT-brain-33c-xbirch-gtf_DGE_lfc-shr_all.csv_with-F1info.csv",header=T)
f1_22c_brain_misexp <- subset(f1_22c_brain, padj_res.22CbirVF1 < 0.1 & padj_res.22CF1Vmal < 0.1 & ((LFC_res.22CbirVF1 > 0 & LFC_res.22CF1Vmal < 0) | (LFC_res.22CbirVF1 < 0 & LFC_res.22CF1Vmal > 0))) #521
f1_33c_brain_misexp <- subset(f1_33c_brain, padj_res.33CbirVF1 < 0.1 & padj_res.33CF1Vmal < 0.1 & ((LFC_res.33CbirVF1 > 0 & LFC_res.33CF1Vmal < 0) | (LFC_res.33CbirVF1 < 0 & LFC_res.33CF1Vmal > 0))) #144
intersect(f1_22c_brain_misexp$Gene,f1_33c_brain_misexp$Gene) #17
f1misexp <- union(f1_22c_brain_misexp$Gene,f1_33c_brain_misexp$Gene)
#521+144-17=648 genes are misexpressed in f1s under either 22c, 33c, or both

# get the intersect of all 33v22c genes in xbir/xmal (134) and the misexxpressed set
intersect(dge33v22_xbir_xmal,f1misexp) #136

# overall, 136/1503 (9.0%) genes that were responsive to temperature in xbirch and/or xmal were misexpressed in f1s under at least one thermal context


### LIVER RNAseq (X. birchmanni pseudoalignment)

### LIVER 33cV22c
xbir_liver_dge_33v22c <- read.csv("input_files/TT-liver-22cV33c-xbirch-gtf_DGE_lfc-shr_all.csv")
xbir_liver_sig <- subset(xbir_liver_dge_33v22c,padj_res.bir33cV22c < 0.1)
f1_liver_sig <- subset(xbir_liver_dge_33v22c,padj_res.F133cV22c < 0.1)
xmal_liver_sig <- subset(xbir_liver_dge_33v22c,padj_res.mal33cV22c < 0.1)
dim(xbir_liver_sig) # 113
dim(f1_liver_sig)   # 408
dim(xmal_liver_sig) # 38

# genes uniquely differentially expressed in response to thermal stress per group
all_sig <- subset(xbir_liver_dge_33v22c,padj_res.bir33cV22c < 0.1 & padj_res.F133cV22c < 0.1 & padj_res.mal33cV22c < 0.1)
bir_only_sig <- subset(xbir_liver_dge_33v22c,padj_res.bir33cV22c < 0.1 & !(padj_res.F133cV22c < 0.1) & !(padj_res.mal33cV22c < 0.1))
f1_only_sig <- subset(xbir_liver_dge_33v22c,!(padj_res.bir33cV22c < 0.1) & padj_res.F133cV22c < 0.1 & !(padj_res.mal33cV22c < 0.1))
mal_only_sig <- subset(xbir_liver_dge_33v22c,!(padj_res.bir33cV22c < 0.1) & !(padj_res.F133cV22c < 0.1) & padj_res.mal33cV22c < 0.1)
bir_mal_only_sig  <- subset(xbir_liver_dge_33v22c,padj_res.bir33cV22c < 0.1 & !(padj_res.F133cV22c < 0.1) & padj_res.mal33cV22c < 0.1)
bir_mal_sig  <- subset(xbir_liver_dge_33v22c,padj_res.bir33cV22c < 0.1 & padj_res.mal33cV22c < 0.1)
bir_mal_sig_negLFC <- subset(xbir_liver_dge_33v22c,padj_res.bir33cV22c < 0.1 & padj_res.mal33cV22c < 0.1 & LFC_res.bir33cV22c < 0 & LFC_res.mal33cV22c < 0)
bir_mal_sig_posLFC <- subset(xbir_liver_dge_33v22c,padj_res.bir33cV22c < 0.1 & padj_res.mal33cV22c < 0.1 & LFC_res.bir33cV22c > 0 & LFC_res.mal33cV22c > 0)
subset(xbir_liver_dge_33v22c,padj_res.bir33cV22c < 0.1 & padj_res.mal33cV22c < 0.1 & LFC_res.bir33cV22c < 0 & LFC_res.mal33cV22c > 0) #0
subset(xbir_liver_dge_33v22c,padj_res.bir33cV22c < 0.1 & padj_res.mal33cV22c < 0.1 & LFC_res.bir33cV22c > 0 & LFC_res.mal33cV22c < 0) #0

dim(all_sig) # 8
dim(bir_only_sig) # 72
dim(f1_only_sig)  # 383
dim(mal_only_sig) # 11
dim(bir_mal_only_sig) # 7
dim(bir_mal_sig) # 17 = all 17 of these are dge in the same direction.
dim(bir_mal_sig_negLFC) # 4
dim(bir_mal_sig_posLFC) # 13


### LIVER 22c
liver_dge_22c_og <- read.csv("Scripts/input_files/TT-liver-22c-xbirch-gtf_DGE_lfc-shr_all.csv_with-F1info.csv")
# first subset by genes that are expressed in both bir and mal under both thermal contexts
liver_dge_22c <- subset(liver_dge_22c_og,!is.na(padj_res.22CbirVmal))
dim(liver_dge_22c) # 17440

## liver 22c: how many genes that are sig diff expressed between xbirVxmal have xbir-like f1 expression?
liver_dge_22c_f1bir <- subset(liver_dge_22c, padj_res.22CbirVmal < 0.1 & !(padj_res.22CbirVF1 < 0.1) & padj_res.22CF1Vmal < 0.1 ) # f1 ge bir-like = 204
dim(subset(liver_dge_22c_f1bir,LFC_res.22CbirVmal < 0 & LFC_res.22CF1Vmal < 0)) # 73 xbir-like f1s have sig lower exp compared to xmal  [f1&bir < mal]
dim(subset(liver_dge_22c_f1bir,LFC_res.22CbirVmal > 0 & LFC_res.22CF1Vmal > 0)) # 131 xbir-like f1s have sig higher exp compared to xmal [f1&bir > mal]
# total: 73+131=204

## liver 22c: how many genes that are sig diff expressed between xbirVxmal have xbir-like f1 expression?
liver_dge_22c_f1mal <- subset(liver_dge_22c, padj_res.22CbirVmal < 0.1 & padj_res.22CbirVF1 < 0.1 & !(padj_res.22CF1Vmal < 0.1) ) # f1 ge mal-like = 135
dim(subset(liver_dge_22c_f1mal,LFC_res.22CbirVmal < 0 & LFC_res.22CbirVF1 < 0)) # 74 xmal-like f1s have sig higher exp compared to xbir [f1&mal > bir]
dim(subset(liver_dge_22c_f1mal,LFC_res.22CbirVmal > 0 & LFC_res.22CbirVF1 > 0)) # 61 xmal-like f1s have sig lower exp compared to xbir  [f1&mal < bir]
# total: 74+61=135

## liver 22c: how many genes that are sig diff expressed for all comparisons have:
liver_dge_22c_allsig <- subset(liver_dge_22c, padj_res.22CbirVmal < 0.1 & padj_res.22CbirVF1 < 0.1 & padj_res.22CF1Vmal < 0.1 ) # 22c all groups sig = 7
#intermediate f1 expression? 7
dim(subset(liver_dge_22c_allsig,LFC_res.22CbirVmal > 0 & LFC_res.22CbirVF1 > 0 & LFC_res.22CF1Vmal > 0 )) # bir > f1 > mal = 4
dim(subset(liver_dge_22c_allsig,LFC_res.22CbirVmal < 0 & LFC_res.22CbirVF1 < 0 & LFC_res.22CF1Vmal < 0 )) # bir < f1 < mal = 3
#sig high f1 misexpression? 0
dim(subset(liver_dge_22c_allsig, LFC_res.22CbirVF1 < 0 & LFC_res.22CF1Vmal > 0 )) # bir&mal < f1 = 0
#sig low f1 misexpression? 0
dim(subset(liver_dge_22c_allsig, LFC_res.22CbirVF1 > 0 & LFC_res.22CF1Vmal < 0 )) # bir&mal > f1 = 0

## liver 22c: how many genes are misexpressed in f1s? 10 / 24173
f1_sig <- subset(liver_dge_22c, padj_res.22CbirVF1 < 0.1 & padj_res.22CF1Vmal < 0.1 )
dim(subset(f1_sig, LFC_res.22CbirVF1 < 0 & LFC_res.22CF1Vmal > 0 )) # bir&mal < f1 = 3
dim(subset(f1_sig, LFC_res.22CbirVF1 > 0 & LFC_res.22CF1Vmal < 0 )) # bir&mal > f1 = 7


### LIVER 33c
liver_dge_33c_og <- read.csv("Scripts/input_files/TT-liver-33c-xbirch-gtf_DGE_lfc-shr_all.csv_with-F1info.csv")
# first subset by genes that are expressed in both bir and mal under both thermal contexts
liver_dge_33c <- subset(liver_dge_33c_og,!is.na(padj_res.33CbirVmal))
dim(liver_dge_33c) # 17002

## liver 33c: how many genes that are sig diff expressed between xbirVxmal have xbir-like f1 expression?
liver_dge_33c_f1bir <- subset(liver_dge_33c, padj_res.33CbirVmal < 0.1 & !(padj_res.33CbirVF1 < 0.1) & padj_res.33CF1Vmal < 0.1 ) # f1 ge bir-like = 79
dim(subset(liver_dge_33c_f1bir,LFC_res.33CbirVmal < 0 & LFC_res.33CF1Vmal < 0)) # 33 xbir-like f1s have sig lower exp compared to xmal  [f1&bir < mal]
dim(subset(liver_dge_33c_f1bir,LFC_res.33CbirVmal > 0 & LFC_res.33CF1Vmal > 0)) # 46 xbir-like f1s have sig higher exp compared to xmal [f1&bir > mal]
# total: 33+46=79

## liver 33c: how many genes that are sig diff expressed between xbirVxmal have xbir-like f1 expression?
liver_dge_33c_f1mal <- subset(liver_dge_33c, padj_res.33CbirVmal < 0.1 & padj_res.33CbirVF1 < 0.1 & !(padj_res.33CF1Vmal < 0.1) ) # f1 ge mal-like = 799
dim(subset(liver_dge_33c_f1mal,LFC_res.33CbirVmal < 0 & LFC_res.33CbirVF1 < 0)) # 410 xmal-like f1s have sig higher exp compared to xbir [f1&mal > bir]
dim(subset(liver_dge_33c_f1mal,LFC_res.33CbirVmal > 0 & LFC_res.33CbirVF1 > 0)) # 389 xmal-like f1s have sig lower exp compared to xbir  [f1&mal < bir]
# total: 410+389=799

## liver 33c: how many genes that are sig diff expressed for all comparisons have:
liver_dge_33c_allsig <- subset(liver_dge_33c, padj_res.33CbirVmal < 0.1 & padj_res.33CbirVF1 < 0.1 & padj_res.33CF1Vmal < 0.1 ) # 33c all groups sig = 29
#intermediate f1 expression? 23
dim(subset(liver_dge_33c_allsig,LFC_res.33CbirVmal > 0 & LFC_res.33CbirVF1 > 0 & LFC_res.33CF1Vmal > 0 )) # bir > f1 > mal = 14
dim(subset(liver_dge_33c_allsig,LFC_res.33CbirVmal < 0 & LFC_res.33CbirVF1 < 0 & LFC_res.33CF1Vmal < 0 )) # bir < f1 < mal = 9
#sig high f1 misexpression? 3
dim(subset(liver_dge_33c_allsig, LFC_res.33CbirVF1 < 0 & LFC_res.33CF1Vmal > 0 )) # bir&mal < f1 = 3
#sig low f1 misexpression? 3
dim(subset(liver_dge_33c_allsig, LFC_res.33CbirVF1 > 0 & LFC_res.33CF1Vmal < 0 )) # bir&mal > f1 = 3

## liver 33c: how many genes are misexpressed in f1s? 96 / 24173
f1_sig <- subset(liver_dge_33c, padj_res.33CbirVF1 < 0.1 & padj_res.33CF1Vmal < 0.1 )
dim(subset(f1_sig, LFC_res.33CbirVF1 < 0 & LFC_res.33CF1Vmal > 0 )) # bir&mal < f1 = 59
dim(subset(f1_sig, LFC_res.33CbirVF1 > 0 & LFC_res.33CF1Vmal < 0 )) # bir&mal > f1 = 37

## liver misexpression summary
# total genes: 17440
# 22C: 10/17440 = 0.0006 = 0.06%
# 33C: 96/17440 = 0.0055 = 0.55%
# 105/17440 = 0.60% genes are misexpressed in f1s under at least one thermal context
# 4/134 = 3.00% genes that were responsive to temperature in mal and/or birch were misexpressed in F1s under at least one thermal context

## how many of the genes with 33cV22c dge in Xbirch and/or Xmal are misexpressed in F1s under at least one thermal context?
# intersect of genes dge between 22v33c brich and 22v33c mal
# liver
exp22v33c_liver <- read.csv("input_files/TT-liver-22cV33c-xbirch-gtf_DGE_lfc-shr_all.csv",header=T)

sig22v33_bir <- subset(exp22v33c_liver, padj_res.bir33cV22c < 0.1) #113
sig22v33_mal <- subset(exp22v33c_liver, padj_res.mal33cV22c < 0.1) #38
intersect(sig22v33_bir$Gene,sig22v33_mal$Gene) #17
dge33v22_xbir_xmal <- union(sig22v33_bir$Gene,sig22v33_mal$Gene)
# 113+38-17=134 genes are dge in birch and/or mal between 22c vs 33c

# pull all misexpressed genes from both thermal contexts
f1_22c_liver_misexp <- subset(liver_dge_22c, padj_res.22CbirVF1 < 0.1 & padj_res.22CF1Vmal < 0.1 & ((LFC_res.22CbirVF1 > 0 & LFC_res.22CF1Vmal < 0) | (LFC_res.22CbirVF1 < 0 & LFC_res.22CF1Vmal > 0))) #10
f1_33c_liver_misexp <- subset(liver_dge_33c, padj_res.33CbirVF1 < 0.1 & padj_res.33CF1Vmal < 0.1 & ((LFC_res.33CbirVF1 > 0 & LFC_res.33CF1Vmal < 0) | (LFC_res.33CbirVF1 < 0 & LFC_res.33CF1Vmal > 0))) #96
intersect(f1_22c_brain_misexp$Gene,f1_33c_brain_misexp$Gene) #1
f1misexp <- union(f1_22c_liver_misexp$Gene,f1_33c_liver_misexp$Gene)
#10+96-1=105 genes are misexpressed in f1s under either 22c, 33c, or both

# get the intersect of all 33v22c genes in xbir/xmal (134) and the misexxpressed set
intersect(dge33v22_xbir_xmal,f1misexp) #4

# overall, 4/134 (3.0%) genes that were responsive to temperature in xbirch and/or xmal were misexpressed in f1s under at least one thermal context

