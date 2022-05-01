## simulate_null_F1_expression.R
##
## Simulate F1 RNAseq transcript counts based on observed parental
## transcript counts.
## 
# MS & CYP I-2022

#install.packages("truncnorm")
library("truncnorm")

## per sample counts to tpm conversion
## takes vector of counts and effective lengths for one sample
## checkout https://support.bioconductor.org/p/91218/
counts2tpm <- function(counts,len) {
  x <- counts/len
  return(x*1e6/sum(x))
}

# use kallisto abundance files to get observed parental expression counts per transcript
# xbir annotations
setwd("input_files/kallisto_output/kallisto_posttrim_xbirch-inferred-txtome")
# 22c brain
mal1<-read.table(file="S1_TSR1-22C-Xmal-brain_1-A1_kallisto/abundance.tsv",head=TRUE)
mal2<-read.table(file="S48_TSR3-22C-Xmal-brain_1-E7_kallisto/abundance.tsv",head=TRUE)
mal3<-read.table(file="S25_TSR2-22C-Xmal-brain_1-F4_kallisto/abundance.tsv",head=TRUE)
bir1<-read.table(file="S5_TSR1-22C-Xbirch-brain_5-E1_kallisto/abundance.tsv",head=TRUE)
bir2<-read.table(file="S29_TSR2-22C-Xbirch-brain_5-B5_kallisto/abundance.tsv",head=TRUE)
bir3<-read.table(file="S52_TSR3-22C-Xbirch-brain_5-A10_kallisto/abundance.tsv",head=TRUE)
# 33c brain
mal1<-read.table(file="S13_TSR1-33C-Xmal-brain_13-E2_kallisto/abundance.tsv",head=TRUE)
mal2<-read.table(file="S36_TSR2-33C-Xmal-brain_13-A6_kallisto/abundance.tsv",head=TRUE)
mal3<-read.table(file="S60_TSR3-33C-Xmal-brain_13-E8_kallisto/abundance.tsv",head=TRUE)
bir1<-read.table(file="S17_TSR1-33C-Xbirch-brain_17-A3_kallisto/abundance.tsv",head=TRUE)
bir2<-read.table(file="S40_TSR2-33C-Xbirch-brain_17-E6_kallisto/abundance.tsv",head=TRUE)
bir3<-read.table(file="S64_TSR3-33C-Xbirch-brain_17-A9_kallisto/abundance.tsv",head=TRUE)

# 22c liver
mal1<-read.table(file="S28_TSR2-22C-Xmal-liver_4-A5_kallisto/abundance.tsv",head=TRUE)
mal2<-read.table(file="S4_TSR1-22C-Xmal-liver_4-B4_kallisto/abundance.tsv",head=TRUE)
mal3<-read.table(file="S51_TSR3-22C-Xmal-liver_4-H7_kallisto/abundance.tsv",head=TRUE)
bir1<-read.table(file="S32_TSR2-22C-Xbirch-liver_8-E5_kallisto/abundance.tsv",head=TRUE)
bir2<-read.table(file="S55_TSR3-22C-Xbirch-liver_8-D10_kallisto/abundance.tsv",head=TRUE)
bir3<-read.table(file="S8_TSR1-22C-Xbirch-liver_8-H1_kallisto/abundance.tsv",head=TRUE)
# 33c liver
mal1<-read.table(file="S16_TSR1-33C-Xmal-liver_16-H2_kallisto/abundance.tsv",head=TRUE)
mal2<-read.table(file="S63_TSR3-33C-Xmal-liver_16-H8_kallisto/abundance.tsv",head=TRUE)
mal3<-read.table(file="S63_TSR3-33C-Xmal-liver_16-H8_kallisto/abundance.tsv",head=TRUE) # no mal3, reuse mal2
bir1<-read.table(file="S20_TSR1-33C-Xbirch-liver_20-D3_kallisto/abundance.tsv",head=TRUE)
bir2<-read.table(file="S43_TSR2-33C-Xbirch-liver_20-H6_kallisto/abundance.tsv",head=TRUE)
bir3<-read.table(file="S67_TSR3-33C-Xbirch-liver_20-D9_kallisto/abundance.tsv",head=TRUE)

# clear the simulated_counts matrices
simulated_counts_f1_1<-{}
simulated_counts_f1_2<-{}
simulated_counts_f1_3<-{}

## simulating estimated counts
# null F1 values generated for each transcript (24475 total)
for(x in 1:length(mal1[,1])){

  # collect observed counts values from parents
  mal_exp<-c(mal1$est_counts[x],mal2$est_counts[x],mal3$est_counts[x])
  bir_exp<-c(bir1$est_counts[x],bir2$est_counts[x],bir3$est_counts[x])

  # get mean and standard deviation of counts per parent
  avg_mal<-mean(mal_exp)
  sd_mal<-sd(mal_exp)

  avg_bir<-mean(bir_exp)
  sd_bir<-sd(bir_exp)

  # get the average length of the transcript over all parent replicates
  eff_length_avg<-mean(c(mal1$eff_length[x],mal2$eff_length[x],mal3$eff_length[x],bir1$eff_length[x],bir2$eff_length[x],bir3$eff_length[x]))

  # for each parent haplotype, create random distribution with center of average count and standard deviation of average count sd
  # ranging from 0 to the maximum counts observed among all parents (bir and mal)
  # pull three values, for three pairs of parents
  ## thought: should there be an upper bound? from all est_counts not for current gene
  mal_hap<-rtruncnorm(3, 0, max(mal1$est_counts,mal2$est_counts,mal3$est_counts,bir1$est_counts,bir2$est_counts,bir3$est_counts), mean = avg_mal, sd = sd_mal)
  bir_hap<-rtruncnorm(3, 0, max(mal1$est_counts,mal2$est_counts,mal3$est_counts,bir1$est_counts,bir2$est_counts,bir3$est_counts), mean = avg_bir, sd = sd_bir)
  mal_hap[is.na(mal_hap)] <- 0
  bir_hap[is.na(bir_hap)] <- 0

  # set F1 count to average of parents
  hyb_val<-rowMeans(cbind(mal_hap,bir_hap))

  # create a matrix of transcript ID, transcript length, and estimated counts per F1
  simulated_counts_f1_1<-rbind(simulated_counts_f1_1,cbind(as.character(mal1$target_id[x]),mal1$length[x],eff_length_avg,hyb_val[1]))
  simulated_counts_f1_2<-rbind(simulated_counts_f1_2,cbind(as.character(mal1$target_id[x]),mal1$length[x],eff_length_avg,hyb_val[2]))
  simulated_counts_f1_3<-rbind(simulated_counts_f1_3,cbind(as.character(mal1$target_id[x]),mal1$length[x],eff_length_avg,hyb_val[3]))

}

simulated_counts_f1_1 <- as.data.frame(simulated_counts_f1_1)
simulated_counts_f1_1$tpm <- counts2tpm(as.numeric(simulated_counts_f1_1$V4),as.numeric(simulated_counts_f1_1$eff_length_avg))
simulated_counts_f1_2 <- as.data.frame(simulated_counts_f1_2)
simulated_counts_f1_2$tpm <- counts2tpm(as.numeric(simulated_counts_f1_2$V4),as.numeric(simulated_counts_f1_2$eff_length_avg))
simulated_counts_f1_3 <- as.data.frame(simulated_counts_f1_3)
simulated_counts_f1_3$tpm <- counts2tpm(as.numeric(simulated_counts_f1_3$V4),as.numeric(simulated_counts_f1_3$eff_length_avg))

colnames(simulated_counts_f1_1)<-colnames(mal1)
colnames(simulated_counts_f1_2)<-colnames(mal1)
colnames(simulated_counts_f1_3)<-colnames(mal1)

# write out
temp <- "33C"
tissue <- "brain"
write.table(simulated_counts_f1_1,paste0("../simulatedF1_kallisto_posttrim_xbirch-inferred-txtome/F1_1_simulated_",temp,"-malxbirchF1-",tissue,"_kallisto/abundance.tsv"),sep="\t",quote=F,row.names=F)
write.table(simulated_counts_f1_2,paste0("../simulatedF1_kallisto_posttrim_xbirch-inferred-txtome/F1_2_simulated_",temp,"-malxbirchF1-",tissue,"_kallisto/abundance.tsv"),sep="\t",quote=F,row.names=F)
write.table(simulated_counts_f1_3,paste0("../simulatedF1_kallisto_posttrim_xbirch-inferred-txtome/F1_3_simulated_",temp,"-malxbirchF1-",tissue,"_kallisto/abundance.tsv"),sep="\t",quote=F,row.names=F)




#################
## NOT USED (but could be useful in the future) 
#################
## simulating TPM
# fill in as desired read counts for hybrids
num_million_reads=c(30e6,35e6,40e6)
# clear the simulated_counts matrices
simulated_counts_f1_1_tpm<-{}
simulated_counts_f1_2_tpm<-{}
simulated_counts_f1_3_tpm<-{}

# null F1 values generated for each transcript (19249 total)
for(x in 1:length(mal1[,1])){

  # collect observed TPM values from parents
  mal_exp<-c(mal1$tpm[x],mal2$tpm[x],mal3$tpm[x])
  bir_exp<-c(bir1$tpm[x],bir2$tpm[x],bir3$tpm[x])

  # get mean and standard deviation of TPM per parent
  avg_mal<-mean(mal_exp)
  sd_mal<-sd(mal_exp)

  avg_bir<-mean(bir_exp)
  sd_bir<-sd(bir_exp)

  # get the average length of the transcript over all parent replicates
  eff_length_avg<-mean(c(mal1$eff_length[x],mal2$eff_length[x],mal3$eff_length[x],bir1$eff_length[x],bir2$eff_length[x],bir3$eff_length[x]))

  # for each parent, create random distribution with center of average tpm and standard deviation of average tpm sd
  # ranging from 0 to the maximum tpm observed among all parents (bir and mal)
  # pull three values, for three pairs of parents
  mal_hap<-rtruncnorm(3, 0, max(mal1$tpm,mal2$tpm,mal3$tpm,bir1$tpm,bir2$tpm,bir3$tpm), mean = avg_mal, sd = sd_mal)
  bir_hap<-rtruncnorm(3, 0, max(mal1$tpm,mal2$tpm,mal3$tpm,bir1$tpm,bir2$tpm,bir3$tpm), mean = avg_bir, sd = sd_bir)
  mal_hap[is.na(mal_hap)] <- 0
  bir_hap[is.na(bir_hap)] <- 0

  # set F1 tpm to average of parents
  hyb_val<-rowMeans(cbind(mal_hap,bir_hap))

  # create a matrix of transcript ID, transcript length, F1 estimated counts, and F1 tpm per F1
  # tpm = (counts / lengths) / sum(counts / lengths) * 1e6
  # counts = tpm
  # }
  simulated_counts_f1_1_tpm<-rbind(simulated_counts_f1_1_tpm,cbind(as.character(mal1$target_id[x]),mal1$length[x],eff_length_avg,(eff_length_avg/1000)*hyb_val[1]*(num_million_reads[1]/1e6),hyb_val[1]))

  simulated_counts_f1_2_tpm<-rbind(simulated_counts_f1_2_tpm,cbind(as.character(mal1$target_id[x]),mal1$length[x],eff_length_avg,(eff_length_avg/1000)*hyb_val[2]*(num_million_reads[2]/1e6),hyb_val[2]))

  simulated_counts_f1_3_tpm<-rbind(simulated_counts_f1_3_tpm,cbind(as.character(mal1$target_id[x]),mal1$length[x],eff_length_avg,(eff_length_avg/1000)*hyb_val[3]*(num_million_reads[3]/1e6),hyb_val[3]))

}

colnames(simulated_counts_f1_1_tpm)<-colnames(mal1)
colnames(simulated_counts_f1_2_tpm)<-colnames(mal1)
colnames(simulated_counts_f1_3_tpm)<-colnames(mal1)



