#load empirical data

data<-read.csv(file="~/Box/Schumer_lab_resources/Project_files/Thermal_tolerance_projects/Data/dge/brain_w-mito/transgressive-exp/TT-brain-33c-w-mito_DGE_lfc-shr_all.csv_with-annots.csv_with-F1info.csv",head=TRUE)
save_genes<-subset(data$Gene,data$F1_transgress.high.57==1)

# testing
#data <- data[1:100,]

#subset F1s, mal, bir for this treatment

F1<-cbind(data$S21,data$S44,data$S68)
xmal<-cbind(data$S13,data$S36,data$S60)
xbir<-cbind(data$S17,data$S40,data$S64)

#set up & run sims

sim_results_add<-{}
sim_results_mald<-{}
sim_results_bird<-{}

# generates:
#   sim_results_add
for(x in 1:length(data[,1])){
sim_mal<-rnorm(3,mean(xmal[x,]),sd(xmal[x,]))
sim_bir<-rnorm(3,mean(xbir[x,]),sd(xbir[x,]))

add_F1<-rowMeans(cbind(sim_mal,sim_bir))
dom_mal<-rnorm(3,mean(xmal[x,]),sd(xmal[x,]))
dom_bir<-rnorm(3,mean(xbir[x,]),sd(xbir[x,]))

testadd<-cbind(c(add_F1,sim_mal,sim_bir),c("f1","f1","f1","mal","mal","mal","bir","bir","bir"))
model<-aov(testadd[,1]~testadd[,2])
sim_results_add<-rbind(sim_results_add,TukeyHSD(model)$testadd[,4])

# Q - should this be dom_mal,sim_mal,sim_mal?
#testmald<-cbind(c(dom_mal,sim_mal,sim_bir),c("f1","f1","f1","mal","mal","mal","bir","bir","bir"))
testmald<-cbind(c(dom_mal,sim_mal,sim_mal),c("f1","f1","f1","mal","mal","mal","bir","bir","bir"))
model<-aov(testmald[,1]~testmald[,2])
sim_results_mald<-rbind(sim_results_mald,TukeyHSD(model)$testmald[,4])

testbird<-cbind(c(dom_bir,sim_bir,sim_bir),c("f1","f1","f1","mal","mal","mal","bir","bir","bir"))
model<-aov(testbird[,1]~testbird[,2])
sim_results_bird<-rbind(sim_results_bird,TukeyHSD(model)$testbird[,4])

}

## write out each results file
rownames(sim_results_add) <- data$Gene
write.csv(sim_results_add, '~/Desktop/sim_results_add.csv')

rownames(sim_results_mald) <- data$Gene
write.csv(sim_results_mald, '~/Desktop/sim_results_mald.csv')

rownames(sim_results_bird) <- data$Gene
write.csv(sim_results_bird, '~/Desktop/sim_results_bird.csv')


##compare nulls to real data

thresh=1e-4
length(subset(data[,1],(data$F1_transgress.high.57==1 | data$F1_transgress.low.26==1) & data$padj_res.33CF1Vmal<thresh & data$padj_res.33CbirVF1<thresh))

length(subset(data[,1],sim_results_add[,1]*length(na.omit(sim_results_add[,1]))<thresh & sim_results_add[,3]*length(na.omit(sim_results_add[,1]))<thresh))

length(subset(data[,1],sim_results_mald[,1]*length(na.omit(sim_results_mald[,1]))<thresh & sim_results_add[,3]*length(na.omit(sim_results_mald[,1]))<thresh))

length(subset(data[,1],sim_results_bird[,1]*length(na.omit(sim_results_bird[,1]))<thresh & sim_results_add[,3]*length(na.omit(sim_results_bird[,1]))<thresh))

