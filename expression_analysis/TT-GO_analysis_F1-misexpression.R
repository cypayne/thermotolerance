## Gene Ontology enrichment analysis: F1 misexpression in the brain
## takes in deseq2 dge output for all genes to build universe
## tests for GO pathway enrichment

library("GOstats")
library("GSEABase")
library("biomaRt")

## TT brain transgressive
c22<-read.csv(file="input_files/TT-brain-22c-w-mito_DGE_lfc-shr_all.csv_with-annots.csv_with-F1info.csv",head=TRUE)
dgeres<-c22
c33<-read.csv(file="input_files/TT-brain-33c-w-mito_DGE_lfc-shr_all.csv_with-annots.csv_with-F1info.csv",head=TRUE)
dgeres<-c33
c33v22<-read.csv(file="input_files/TT-brain-33c-w-mito_DGE_lfc-shr_all.csv_with-annots.csv_with-F1info.csv",head=TRUE)

mart <- useMart(biomart = "ensembl", dataset = "xmaculatus_gene_ensembl", host="uswest.ensembl.org")
#attributes <- listAttributes(mart)
#attributes[1:50,]

# match go ids to ensembl gene ids
results <- getBM(attributes = c("go_id","external_gene_name","ensembl_gene_id","kegg_enzyme"), filters=c("ensembl_gene_id"),values=dgeres$Gene, mart = mart)

# subset gene universe to only include genes with valid go ids and external gene names
gene_universe<-subset(results,nchar(results$go_id)>0 & nchar(results$external_gene_name) > 0)
gene_universe$ensembl_id<-gene_universe[,3]
gene_universe[,3]<-as.numeric(as.factor(gene_universe[,3]))
gene_universe$Evidence<-rep("ISA",length(gene_universe[,3]))
colnames(gene_universe)<-c("frame.go_id","frame.gene_name","frame.gene_id","frame.KEGG","frame.ensembl","frame.Evidence")

goframeData <- data.frame(gene_universe$frame.go_id,gene_universe$frame.Evidence,gene_universe$frame.gene_id)

goFrame <- GOFrame(goframeData,organism="Xiphophorus")
goAllFrame <- GOAllFrame(goFrame)
gsc <- GeneSetCollection(goAllFrame, setType = GOCollection())

universe <- goframeData$gene_universe.frame.gene_id

## GO enrichment

## subset: TT brain misexpressed f1 genes
# to get immune genes not in misexp set: enrich for all !(misexpressed)
notmisexp.22c <- subset(c22, is.na(sig.F1_exp_profile))
dgesig<-notmisexp.22c

# 22c
head(dgeres)
trans.high.22c<-subset(c22, F1_transgress.high.111 == 1)
dgesig<-trans.high.22c
# trans.low<-subset(dgeres, F1_transgress.low.114 == 1)
# dgesig<-trans.low

# 33c
head(dgeres)
trans.high.33c<-subset(c33, F1_transgress.high.57 == 1)
#dgesig <- trans.high.33c
#trans.low<-subset(dgeres, F1_transgress.low.26 == 1)
#dgesig<-trans.low

# trans.high at both 22c and 33c
both.high <- trans.high.22c[trans.high.22c$Gene %in% trans.high.33c$Gene,]
both.high # just 2 conserved between conditions,
# sqrdl/SQOR: sulfide quinone oxidoreductase, decreases toxic conc. of sulfide [mitochondria]
# hgd: enzyme homogentisate 1,2 dioxygenase, catabolizes tyrosine and phenylalanine [metabolism]

#read in significant genes and match them to the gene_universe ids
#dgesig<-read.csv(file="/Users/cypayne/Desktop/Schumer_lab/sword-qtl/sword-regen-paper/ALd10-2birtrim-allchr_DGE_lfc-shr_res.d10malVbir.padj0.1.csv",head=TRUE)
#dgesig<-read.csv(file="/Users/cypayne/Desktop/Schumer_lab/sword-qtl/sword-regen-paper/ALd10-2maltrim-allchr_DGE_lfc-shr_res.d10malVbir.padj0.1.csv",head=TRUE)

# ## TT: subset padj < 0.1
# # brain bir 33c v 22c
# dim(bir) # 19106
# dgesig <- subset(bir, padj < 0.1)
# dim(dgesig) # 759
# 
# # brain mal 33c v 22c
# dim(mal) # 19106
# dgesig <- subset(mal, padj < 0.1)
# dim(dgesig) # 775
# 
# # liver bir 33c v 22c
# dim(bir) # 19106
# dgesig <- subset(bir, padj < 0.1)
# dim(dgesig) # 72
# 
# # liver mal 33c v 22c
# dim(mal) # 19106
# dgesig <- subset(mal, padj < 0.1)
# dim(dgesig) # 27


genes_sig <- dgesig$Gene
genes_match<-gene_universe[gene_universe$frame.ensembl %in% genes_sig,]
genes_match_sig <- genes_match$frame.gene_id

c22low_go2gene_matchup <- genes_match
write.csv(c22low_go2gene_matchup,'TT-brain-22c-w-mito_DGE_lfc-shr_all-non-misexp-GOsiggenesmatch.csv')

# set params for hyperGTest
# testDirection = "over" for overrepresented genes, = "under" for underrepresented
# interested in overrep genes
# three categories: 
#   cellular component (CC; where gene products are active)
#   molecular function (MF; the biological function of gene or gene product) 
#   biological process (BP; pathways or larger processes that multiple gene products involved in).
# Output: 
# ExpCount is the expected count and the Count is how many instances of that term were actually oberved 
# in your gene list while the Size is the number that could have been found in your gene list if every 
# instance had turned up. Values like the ExpCount and the Size are going to be affected by what is included 
# in the gene universe as well as by whether or not it was a conditional test.
params <- GSEAGOHyperGParams(name="Xiphophorus maculatus genes",
                             geneSetCollection=gsc,
                             geneIds = genes_match_sig,
                             universeGeneIds = universe,
                             ontology = "BP",
                             pvalueCutoff = 0.95,
                             conditional = FALSE,
                             testDirection = "over")

OverBP <- hyperGTest(params)
results_OverBP<-summary(OverBP)

## Output an additional column with the genes falling under each GO category
# FOR BP
res_GOBPID <- results_OverBP$GOBPID # for BP
cats <- geneIdsByCategory(OverBP) 
# get list of lists: all geneids per GOBPID
extract_geneids <- sapply(res_GOBPID, function(go) cats[[go]])
# this gives dataframe of two columns: GOBPID, and the list of gene ids falling under it
gobpid2geneid <- data.frame(GOBPID=names(extract_geneids), GO_GENEID_list=matrix(extract_geneids))
gobpid2geneid[sapply(gobpid2geneid, is.list)] <- apply(gobpid2geneid[sapply(gobpid2geneid, is.list)], 1, function(x) paste(unlist(x), sep=", ", collapse=", "))
# merge GO_GENEID_list with the rest of GO results
go_output <- merge(results_OverBP,gobpid2geneid, by.x="GOBPID", by.y="GOBPID", all=TRUE)
write.table(go_output,'TT-brain-22c-w-mito_DGE_lfc-shr_all-non-misexp-GOresults-BP_w-geneids-pval0.95.tsv', row.names = FALSE, quote=FALSE, sep="\t")

# FOR CC
res_GOCCID <- results_OverBP$GOCCID # for CC
cats <- geneIdsByCategory(OverBP) 
# get list of lists: all geneids per GOBPID
extract_geneids <- sapply(res_GOCCID, function(go) cats[[go]])
# this gives dataframe of two columns: GOBPID, and the list of gene ids falling under it
goccid2geneid <- data.frame(GOCCID=names(extract_geneids), GO_GENEID_list=matrix(extract_geneids))
goccid2geneid[sapply(goccid2geneid, is.list)] <- apply(goccid2geneid[sapply(goccid2geneid, is.list)], 1, function(x) paste(unlist(x), sep=", ", collapse=", "))
# merge GO_GENEID_list with the rest of GO results
go_output <- merge(results_OverBP,goccid2geneid, by.x="GOCCID", by.y="GOCCID", all=TRUE)
write.table(go_output,'TT-brain-22c-w-mito_DGE_lfc-shr_all-F1trans.high-GOresults-CC_w-geneids.tsv', row.names = FALSE, quote=FALSE, sep="\t")

# 33c
#write.csv(results_OverBP,'/Users/cypayne/Box/Schumer_lab_resources/Project_files/Thermal_tolerance_projects/Data/dge/brain_w-mito/transgressive-exp/TT-brain-33c-w-mito_DGE_lfc-shr_all-F1trans.high-GOresults-BP.csv')
#write.csv(results_OverBP,'/Users/cypayne/Box/Schumer_lab_resources/Project_files/Thermal_tolerance_projects/Data/dge/brain_w-mito/transgressive-exp/TT-brain-33c-w-mito_DGE_lfc-shr_all-F1trans.high-GOresults-CC.csv')
#write.csv(results_OverBP,'/Users/cypayne/Box/Schumer_lab_resources/Project_files/Thermal_tolerance_projects/Data/dge/brain_w-mito/transgressive-exp/TT-brain-33c-w-mito_DGE_lfc-shr_all-F1trans.high-GOresults-MF.csv')
#write.csv(results_OverBP,'/Users/cypayne/Box/Schumer_lab_resources/Project_files/Thermal_tolerance_projects/Data/dge/brain_w-mito/transgressive-exp/TT-brain-33c-w-mito_DGE_lfc-shr_all-F1trans.low-GOresults-BP.csv')
#write.csv(results_OverBP,'/Users/cypayne/Box/Schumer_lab_resources/Project_files/Thermal_tolerance_projects/Data/dge/brain_w-mito/transgressive-exp/TT-brain-33c-w-mito_DGE_lfc-shr_all-F1trans.low-GOresults-CC.csv')
#write.csv(results_OverBP,'/Users/cypayne/Box/Schumer_lab_resources/Project_files/Thermal_tolerance_projects/Data/dge/brain_w-mito/transgressive-exp/TT-brain-33c-w-mito_DGE_lfc-shr_all-F1trans.low-GOresults-MF.csv')


# # 22c
#write.csv(results_OverBP,'/Users/cypayne/Box/Schumer_lab_resources/Project_files/Thermal_tolerance_projects/Data/dge/brain_w-mito/transgressive-exp/TT-brain-22c-w-mito_DGE_lfc-shr_all-F1trans.high-GOresults-BP-2.csv')
# #write.csv(results_OverBP,'/Users/cypayne/Box/Schumer_lab_resources/Project_files/Thermal_tolerance_projects/Data/dge/brain_w-mito/transgressive-exp/TT-brain-22c-w-mito_DGE_lfc-shr_all-F1trans.high-GOresults-CC.csv')
# write.csv(results_OverBP,'/Users/cypayne/Box/Schumer_lab_resources/Project_files/Thermal_tolerance_projects/Data/dge/brain_w-mito/transgressive-exp/TT-brain-22c-w-mito_DGE_lfc-shr_all-F1trans.high-GOresults-MF.csv')
# #write.csv(results_OverBP,'/Users/cypayne/Box/Schumer_lab_resources/Project_files/Thermal_tolerance_projects/Data/dge/brain_w-mito/transgressive-exp/TT-brain-22c-w-mito_DGE_lfc-shr_all-F1trans.low-GOresults-BP.csv')
# #write.csv(results_OverBP,'/Users/cypayne/Box/Schumer_lab_resources/Project_files/Thermal_tolerance_projects/Data/dge/brain_w-mito/transgressive-exp/TT-brain-22c-w-mito_DGE_lfc-shr_all-F1trans.low-GOresults-CC.csv')
# write.csv(results_OverBP,'/Users/cypayne/Box/Schumer_lab_resources/Project_files/Thermal_tolerance_projects/Data/dge/brain_w-mito/transgressive-exp/TT-brain-22c-w-mito_DGE_lfc-shr_all-F1trans.low-GOresults-MF.csv')


# match sig go ids to gene nums
#go2gene_matchup <- data.frame("GOBPID" = genes_match$frame.go_id, "GeneID" = genes_match$frame.ensembl)

#c33high_go2gene_matchup <- genes_match
#write.csv(c33high_go2gene_matchup,'/Users/cypayne/Box/Schumer_lab_resources/Project_files/Thermal_tolerance_projects/Data/dge/brain_w-mito/transgressive-exp/TT-brain-33c-w-mito_DGE_lfc-shr_all-F1trans.high-GOsiggenesmatch.csv')
# c33low_go2gene_matchup <- genes_match
# write.csv(c33low_go2gene_matchup,'/Users/cypayne/Box/Schumer_lab_resources/Project_files/Thermal_tolerance_projects/Data/dge/brain_w-mito/transgressive-exp/TT-brain-33c-w-mito_DGE_lfc-shr_all-F1trans.low-GOsiggenesmatch.csv')

#c22high_go2gene_matchup <- genes_match
#write.csv(c22high_go2gene_matchup,'/Users/cypayne/Box/Schumer_lab_resources/Project_files/Thermal_tolerance_projects/Data/dge/brain_w-mito/transgressive-exp/TT-brain-22c-w-mito_DGE_lfc-shr_all-F1trans.high-GOsiggenesmatch.csv')
c22low_go2gene_matchup <- genes_match
write.csv(c22low_go2gene_matchup,'TT-brain-22c-w-mito_DGE_lfc-shr_all-F1trans.low-GOsiggenesmatch.csv')

# to look at genes that fall under a given go category
# in sig genes
# 48, 49, 108, 126, 130, 139
cats <- geneIdsByCategory(OverBP)
c <- cats$`GO:0000422` # autophagy of mitochondrion:  atg7, rb1cc1
c <- cats$`GO:0061726` # mitochondrion disassembly:   atg7, rb1cc1
c <- cats$`GO:0006839` # mitochondrial transport:     SLC25A30, slc25a38b, chchd4a
c <- cats$`GO:0045041` # protein import into mitochondrial intermembrane space: chchd4a
c <- cats$`GO:0000423` # mitophagy: atg7
c <- cats$`GO:0070125` # mitochondrial translational elongation: gfm1
genes_match[genes_match$frame.gene_id %in% c,]

#library(plyr)
#gobpid2geneid <- plyr::ldply(extract_geneids, cbind) # gives two column df, one-to-one GO:geneid map
# now add xmac gene id to this df

# Look at dif in GO enrichment between pseudoalign refs
mal <- read.csv("TT-brain-mal33cV22c-padj0.1-GOresults.csv", header=TRUE)
bir <- read.csv("TT-brain-bir33cV22c-padj0.1-GOresults.csv", header=TRUE)

intersection <- intersect(bir$Term,mal$Term)
dim(bir) #161
dim(mal) #154
length(intersection) #84

bir_overlap <- subset(bir,Term %in% intersection)
bir_nolap <- subset(bir,!(Term %in% intersection)) #77

mal_overlap <- subset(mal,Term %in% intersection)
mal_nolap <- subset(mal,!(Term %in% intersection)) #70
