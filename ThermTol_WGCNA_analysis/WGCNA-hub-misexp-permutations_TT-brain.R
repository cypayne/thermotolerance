## WGCNA-hub-misexp-permutation_TT-liver.R
##
## This script is broken up into the following functions:
## Permute the expected number of hub genes that are transgressive/misexpressed in F1s
## Permute the expected number of transgressive/misexpressed genes in F1s per module
## Permute the expected number of hub genes per module 
##
## Contains code adapted from 
## Morgan et al 2020: https://doi.org/10.1093/molbev/msaa002
##
## cyp I-2022

library(dplyr)

## Check overlap between kTotal and kWithin-filtered hub gene sets
hub.genes_kTotal<-read.csv('input_files/TT-brain-WGCNA-onehot_hubgenes_kTotal0.95_MM0.85.csv',header=T,sep=',') # using 0.95q kTotal: 705 hub genes 
hub.genes_kWithin<-read.csv('input_files/TT-brain-WGCNA-onehot_hubgenes_kWithin0.95_MM0.85.csv',header=T,sep=',') # using 0.95q kWithin: 600 hub genes 
nrow(subset(hub.genes_kTotal, X %in% hub.genes_kWithin$X)) # 517 genes overlap between these hub gene sets

# use kWithin-filtered hub genes for further analysis
hub.genes<-hub.genes_kWithin$X
num_hub_genes <- length(hub.genes) # 600

## First determine how many hub genes are transgressive (low, high, both) at 22C
## 15
allgenes22C<-read.csv("input_files/TT-brain-22c-w-mito_DGE_lfc-shr_all.csv_with-annots.csv_with-F1info.csv",h=T,sep=",")
trans.both<-subset(allgenes22C,subset=F1_transgress.high=="1"|F1_transgress.low=="1") #225
trans.both_hubgenes<-subset(subset(allgenes22C,subset=F1_transgress.high=="1"|F1_transgress.low=="1"), Gene %in% hub.genes) #15
trans.low_hubgenes<-subset(subset(allgenes22C,subset=F1_transgress.low=="1"), Gene %in% hub.genes) #12
trans.high_hubgenes<-subset(subset(allgenes22C,subset=F1_transgress.high=="1"), Gene %in% hub.genes) #3

## Determine how many hub genes are transgressive (low, high, both) at 33C
## NONE
allgenes33C<-read.csv("input_files/TT-brain-33c-w-mito_DGE_lfc-shr_all.csv_with-annots.csv_with-F1info.csv",h=T,sep=",")
trans.both<-subset(allgenes33C,subset=F1_transgress.high=="1"|F1_transgress.low=="1") #83
trans.both_hubgenes<-subset(subset(allgenes33C,subset=F1_transgress.high=="1"|F1_transgress.low=="1"), Gene %in% hub.genes) #0
trans.low_hubgenes<-subset(subset(allgenes33C,subset=F1_transgress.low=="1"), Gene %in% hub.genes) #0
trans.high_hubgenes<-subset(subset(allgenes33C,subset=F1_transgress.high=="1"), Gene %in% hub.genes) #0


## How many hub genes are expected to have transgressive f1 expression (both low and high)?
## Use permutations to determine what the expected distribution of transgressive hub genes is at 22C
allGenes_withinTransgressive<-c()
for (i in 1:1000) {
  # randomly sample num_hub_genes rows from the GeneSet dataframe 
  random_genes<-sample_n(allgenes22C, num_hub_genes)
  ran<-subset(random_genes,subset=F1_transgress.high=="1"|F1_transgress.low=="1")
  number_transgressive<-nrow(ran)
  allGenes_withinTransgressive<-append(allGenes_withinTransgressive,number_transgressive)
}
hist(allGenes_withinTransgressive)
quantile(allGenes_withinTransgressive,probs=c(0.025,0.975))
#Results: 22C:Kwithin=600 95%CI= 2-13. 15 (the number of hub genes that are transgressive) is outside the CI, i.e. we observe significantly more transgressive hub genes than expected by chance at 22C.

# How many hub genes are expected to have low transgressive f1 expression?
allGenes_withLowTransExp<-c()
for (i in 1:1000) {
  #randomly sample num_hub_genes rows from the GeneSet dataframe
  random_genes<-sample_n(allgenes22C, num_hub_genes)
  ran<-subset(random_genes,subset=F1_transgress.low=="1")
  number_transgressive<-nrow(ran)
  allGenes_withLowTransExp<-append(allGenes_withLowTransExp,number_transgressive)
}
hist(allGenes_withLowTransExp)
quantile(allGenes_withLowTransExp,probs=c(0.025,0.975))
#Results: 95%CI= 1-7. 12 is outside the CI, more hub genes are lowly transgressive than expected by chance.

# How many hub genes are expected to have high transgressive f1 expression?
allGenes_withHighTransExp<-c()
for (i in 1:1000) {
  #randomly sample num_hub_genes rows from the GeneSet dataframe
  random_genes<-sample_n(allgenes22C, num_hub_genes)
  ran<-subset(random_genes,subset=F1_transgress.high=="1")
  number_transgressive<-nrow(ran)
  allGenes_withHighTransExp<-append(allGenes_withHighTransExp,number_transgressive)
}
hist(allGenes_withHighTransExp)
quantile(allGenes_withHighTransExp,probs=c(0.025,0.975))
#Results: 95%CI = 1-8. 3 falls within the CI, the observed number of highly transgressive hub genes can be due to chance.


## How many hub genes are expected to have transgressive f1 expression (both low and high) at 33C?
allGenes_withinTransgressive<-c()
for (i in 1:1000) {
  #randomly sample num_hub_genes rows from the GeneSet dataframe 
  random_genes<-sample_n(allgenes33C, num_hub_genes)
  ran<-subset(random_genes,subset=F1_transgress.high=="1"|F1_transgress.low=="1")
  number_transgressive<-nrow(ran)
  allGenes_withinTransgressive<-append(allGenes_withinTransgressive,number_transgressive)
}
hist(allGenes_withinTransgressive)
quantile(allGenes_withinTransgressive,probs=c(0.025,0.975))
#Results: 22C:Kwithin=600 95%CI = 0-6. 0 (the number of hub genes that are transgressive at 33c) falls within the CI.


### WGCNA module permutations
# load all of the significant module names (sigMEs object)
load(file = "input_files/TT-brain-WGCNA_sigMEs.RData")
sigMEs

### Permutations: How many transgressive genes do we expect in each module?
header="TT-brain-WGCNA-onehot_"
output_file=paste0(header,"num-trans-genes-in-ME_permutation-results.csv")
header_line = data.frame('module_name','num_MEgenes','num_MEtrans.both.22c','num_MEtrans.low.22c','num_MEtrans.high.22c','num_MEtrans.both.33c','num_MEtrans.low.33c','num_MEtrans.high.33c','22c_q0.025','22c_q0.975','33c_q0.025','33c_q0.975')
write.table(header_line,output_file,col.names=F,row.names=F,sep=",")
for( MEname in sigMEs ) { 
  MEname <- substring(MEname, 3)
  MEfile <- read.csv(paste0(header,MEname,"_genes.csv"),header=T)
  genes <- MEfile$x
  num_MEgenes <- length(genes)
  # number of transgressive genes in module
  MEtrans.both22c<-subset(subset(allgenes22C,subset=F1_transgress.high=="1"|F1_transgress.low=="1"), Gene %in% genes)
  MEtrans.low22c<-subset(subset(allgenes22C,subset=F1_transgress.low=="1"), Gene %in% genes)
  MEtrans.high22c<-subset(subset(allgenes22C,subset=F1_transgress.high=="1"), Gene %in% genes)
  
  MEtrans.both33c<-subset(subset(allgenes33C,subset=F1_transgress.high=="1"|F1_transgress.low=="1"), Gene %in% genes)
  MEtrans.low33c<-subset(subset(allgenes33C,subset=F1_transgress.low=="1"), Gene %in% genes)
  MEtrans.high33c<-subset(subset(allgenes33C,subset=F1_transgress.high=="1"), Gene %in% genes)
  
  row2out <- data.frame(MEname,num_MEgenes,dim(MEtrans.both22c)[1],dim(MEtrans.low22c)[1],dim(MEtrans.high22c)[1],dim(MEtrans.both33c)[1],dim(MEtrans.low33c)[1],dim(MEtrans.high33c)[1])
  
  # output the expression profiles for the transgressive genes
  write.csv(MEtrans.both22c,paste0(header,MEname,"_transF1-expression-22c.csv"))
  write.csv(MEtrans.both33c,paste0(header,MEname,"_transF1-expression-33c.csv"))
  
  # perms for 22c
  allGenes_withinTransgressive22c<-c()
  for (i in 1:1000) {
    #randomly sample the number of genes in module 
    random_genes<-sample_n(allgenes22C, num_MEgenes)
    ran<-subset(random_genes,subset=F1_transgress.high=="1"|F1_transgress.low=="1")
    number_transgressive<-nrow(ran)
    allGenes_withinTransgressive22c<-append(allGenes_withinTransgressive22c,number_transgressive)
  }
  #hist(allGenes_withinTransgressive22c)
  quantiles22c<-quantile(allGenes_withinTransgressive22c,probs=c(0.025,0.975))
  row2out<-append(row2out,quantiles22c)
  
  # perms for 33c
  allGenes_withinTransgressive33c<-c()
  for (i in 1:1000) {
    #randomly sample the number of genes in module 
    random_genes<-sample_n(allgenes33C, num_MEgenes)
    ran<-subset(random_genes,subset=F1_transgress.high=="1"|F1_transgress.low=="1")
    number_transgressive<-nrow(ran)
    allGenes_withinTransgressive33c<-append(allGenes_withinTransgressive33c,number_transgressive)
  }
  #hist(allGenes_withinTransgressive33c)
  quantiles33c<-quantile(allGenes_withinTransgressive33c,probs=c(0.025,0.975))
  row2out<-append(row2out,quantiles33c)
  
  # output results for each module to the output file
  write.table(row2out,output_file,col.names=F,row.names=F,sep=",",append=T)
  
}

## Permutations: How many hub genes do we expect in each module?
# --> are there more hub (highly connected) genes in the temperature-associated modules?
# observed value: the number of hub genes in each module
# simulations:
#   randomly sample genes from all genes (same number as in module)
#   see how many of these randomly sampled genes are hub genes

hub.filter <- "kTotal" ## Set filter to use (kTotal or kWithin)
if(hub.filter=="kWithin"){ 
  hub.genes<-hub.genes_kWithin$X #600 
} else if(hub.filter=="kTotal"){ 
  hub.genes<-hub.genes_kTotal$X #705 
}
header="TT-brain-WGCNA-onehot_"
output_file=paste0(header,"num",hub.filter,"-hub-genes-in-ME_permutation-results.csv")
header_line = data.frame('module_name','num_MEgenes','num_MEhub.genes','q0.025','q0.0975')
write.table(header_line,output_file,col.names=F,row.names=F,sep=",")
for( MEname in sigMEs ) { 
  MEname <- substring(MEname, 3)
  MEfile <- read.csv(paste0(header,MEname,"_genes.csv"),header=T)
  MEgenes <- MEfile$x
  num_MEgenes <- length(MEgenes)
  # number of hub genes in module
  MEhub.genes<-subset(subset(allgenes22C, Gene %in% hub.genes), Gene %in% MEgenes)
  row2out <- data.frame(MEname,num_MEgenes,dim(MEhub.genes)[1])
  
  # sample from either 22c or 33c, doesn't change hub gene status
  allGenes_withinHub<-c()
  for (i in 1:1000) {
    #randomly sample the number of genes in module
    random_genes<-sample_n(allgenes22C, num_MEgenes)
    ran<-subset(random_genes,Gene %in% hub.genes)
    number_hub<-nrow(ran)
    allGenes_withinHub<-append(allGenes_withinHub,number_hub)
  }
  quantiles<-quantile(allGenes_withinHub,probs=c(0.025,0.975))
  row2out<-append(row2out,quantiles)
  
  # output results for each module to the output file
  write.table(row2out,output_file,col.names=F,row.names=F,sep=",",append=T)
  
}


## Permutations: How many transgressive hub genes do we expect, overall?
## Run for both kWithin and kTotal-filtered hub genes

hub.filter <- "kTotal" ## Set filter to use (kTotal or kWithin)
if(hub.filter=="kWithin"){ 
  hub.genes<-hub.genes_kWithin$X #600 
} else if(hub.filter=="kTotal"){ 
  hub.genes<-hub.genes_kTotal$X #705 
}

header="TT-brain-WGCNA-onehot_"
output_file=paste0(header,"num-trans-",hub.filter,"-hub-genes_permutation-results.csv")
header_line = data.frame('num_hub_genes','num_hub_trans.both.22c','num_hub_trans.low.22c','num_hub_trans.high.22c','num_hub_trans.both.33c','num_hub_trans.low.33c','num_hub_trans.high.33c','22c_q0.025','22c_q0.975','33c_q0.025','33c_q0.975')
write.table(header_line,output_file,col.names=F,row.names=F,sep=",")

# number of hub genes with transgression in F1s
hub.trans.both22c<-subset(subset(allgenes22C,subset=F1_transgress.high=="1"|F1_transgress.low=="1"), Gene %in% hub.genes)
hub.trans.low22c<-subset(subset(allgenes22C,subset=F1_transgress.low=="1"), Gene %in% hub.genes)
hub.trans.high22c<-subset(subset(allgenes22C,subset=F1_transgress.high=="1"), Gene %in% hub.genes)
  
hub.trans.both33c<-subset(subset(allgenes33C,subset=F1_transgress.high=="1"|F1_transgress.low=="1"), Gene %in% hub.genes)
hub.trans.low33c<-subset(subset(allgenes33C,subset=F1_transgress.low=="1"), Gene %in% hub.genes)
hub.trans.high33c<-subset(subset(allgenes33C,subset=F1_transgress.high=="1"), Gene %in% hub.genes)

row2out <- data.frame(length(hub.genes),dim(hub.trans.both22c)[1],dim(hub.trans.low22c)[1],dim(hub.trans.high22c)[1],dim(hub.trans.both33c)[1],dim(hub.trans.low33c)[1],dim(hub.trans.high33c)[1])
  
# output the expression profiles for the hub genes
write.csv(hub.trans.both22c,paste0(header,"_all-trans-",hub.filter,"-hub-expression-22c.csv"))
write.csv(hub.trans.both33c,paste0(header,"_all-trans-",hub.filter,"-hub-expression-33c.csv"))

# perms for 22c
allGenes_withinTransgressive22c<-c()
for (i in 1:1000) {
  #randomly sample the number of genes in module 
  random_genes<-sample_n(allgenes22C, num_hub_genes)
  ran<-subset(random_genes,subset=F1_transgress.high=="1"|F1_transgress.low=="1")
  number_transgressive<-nrow(ran)
  allGenes_withinTransgressive22c<-append(allGenes_withinTransgressive22c,number_transgressive)
}
#hist(allGenes_withinTransgressive22c)
quantiles22c<-quantile(allGenes_withinTransgressive22c,probs=c(0.025,0.975))
row2out<-append(row2out,quantiles22c)
  
# perms for 33c
allGenes_withinTransgressive33c<-c()
for (i in 1:1000) {
  #randomly sample the number of genes in module 
  random_genes<-sample_n(allgenes33C, num_hub_genes)
  ran<-subset(random_genes,subset=F1_transgress.high=="1"|F1_transgress.low=="1")
  number_transgressive<-nrow(ran)
  allGenes_withinTransgressive33c<-append(allGenes_withinTransgressive33c,number_transgressive)
}
#hist(allGenes_withinTransgressive33c)
quantiles33c<-quantile(allGenes_withinTransgressive33c,probs=c(0.025,0.975))
row2out<-append(row2out,quantiles33c)
  
# output results for each module to the output file
write.table(row2out,output_file,col.names=F,row.names=F,sep=",",append=T)
  
