## Gene Ontology enrichment analysis: TT 22c, 33c, and 33cV22c
## takes in deseq2 dge output for all genes to build universe
##                        and for sig genes to test enrichment
# cyp & ms III-2020

library("GOstats")
library("GSEABase")
library("biomaRt")

## TT brain
brain_22c<-read.csv(file="Scripts/input_files/TT-brain-22c-w-mito_DGE_lfc-shr_all.csv",head=TRUE)
brain_33c<-read.csv(file="Scripts/input_files/TT-brain-33c-w-mito_DGE_lfc-shr_all.csv",head=TRUE)
brain_33cV22c<-read.csv(file="Scripts/input_files/TT-brain-33cV22c-w-mito_DGE_lfc-shr_all.csv",head=TRUE)

## TT liver
liver_22c<-read.csv(file="input_files/TT-liver-22c-w-mito_DGE_lfc-shr_all.csv",head=TRUE)
liver_33c<-read.csv(file="input_files/TT-liver-33c-w-mito_DGE_lfc-shr_all.csv",head=TRUE)

# set the data you want to get GO enrichment for
dgeres <- brain_33cV22c

mart <- useMart(biomart = "ensembl", dataset = "xmaculatus_gene_ensembl", host="uswest.ensembl.org")
#attributes <- listAttributes(mart)
#attributes[1:50,]

# match go ids to ensembl gene ids
results <- getBM(attributes = c("go_id","external_gene_name","ensembl_gene_id"), filters=c("ensembl_gene_id"),values=dgeres$Gene, mart = mart)

# subset gene universe to only include genes with valid go ids and external gene names
gene_universe<-subset(results,nchar(results$go_id)>0 & nchar(results$external_gene_name) > 0)
gene_universe$ensembl_id<-gene_universe[,3]
gene_universe[,3]<-as.numeric(as.factor(gene_universe[,3]))
gene_universe$Evidence<-rep("ISA",length(gene_universe[,3]))
colnames(gene_universe)<-c("frame.go_id","frame.gene_name","frame.gene_id","frame.ensembl","frame.Evidence")

goframeData <- data.frame(gene_universe$frame.go_id,gene_universe$frame.Evidence,gene_universe$frame.gene_id)

goFrame <- GOFrame(goframeData,organism="Xiphophorus")
goAllFrame <- GOAllFrame(goFrame)
gsc <- GeneSetCollection(goAllFrame, setType = GOCollection())

universe <- goframeData$gene_universe.frame.gene_id

## GO enrichment

## brain
sig_thresh <- 0.1
# if 22c birVmal
dgesig <- subset(dgeres,padj_res.22CbirVmal < sig_thresh) # 0.1 = 2536
# if 33c birVmal
dgesig <- subset(dgeres,padj_res.33CbirVmal < sig_thresh) # 0.1 = 2422
# if 33cV22c F1
dgesig <- subset(dgeres, padj_res.F133cV22c < sig_thresh) # 0.1 = 1833
# if 33cV22c bir
dgesig <- subset(dgeres, padj_res.bir33cV22c < sig_thresh) # 0.1 = 758, 0.05 = 515
# if 33cV22c mal
dgesig <- subset(dgeres, padj_res.mal33cV22c < sig_thresh) # 0.1 = 763, 0.05 = 541
dim(dgesig)

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

# optional: write the sig gene GO info to an outfile
#gene2go_matchup_22c_brain <- genes_match
#write.csv(gene2go_matchup_22c_brain,'./TT-brain-22c-w-mito_DGE_lfc-shr_padj0.1-birVmal_gene2GO-matchup.csv')
#gene2go_matchup_33c_brain <- genes_match
#write.csv(gene2go_matchup_33c_brain,'./TT-brain-33c-w-mito_DGE_lfc-shr_padj0.1-birVmal_gene2GO-matchup.csv')
#gene2go_matchup_33cV22c_brain <- genes_match
#write.csv(gene2go_matchup_33cV22c_brain,'./TT-brain-33cV22c-w-mito_DGE_lfc-shr_padj0.1-33v22F1_gene2GO-matchup.csv')
#gene2go_matchup_33cV22c_brain <- genes_match
#write.csv(gene2go_matchup_33cV22c_brain,'./TT-brain-33cV22c-w-mito_DGE_lfc-shr_padj0.1-33v22Bir_gene2GO-matchup.csv')
gene2go_matchup_33cV22c_brain <- genes_match
write.csv(gene2go_matchup_33cV22c_brain,'./TT-brain-33cV22c-w-mito_DGE_lfc-shr_padj0.1-33v22Mal_gene2GO-matchup.csv')

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
                              pvalueCutoff = 0.05,
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
#write.table(go_output,'./TT-brain-22c-w-mito_DGE_lfc-shr_pval0.05_birVmal_BP.tsv', row.names = FALSE, quote=FALSE, sep="\t")
#write.table(go_output,'./TT-brain-33c-w-mito_DGE_lfc-shr_pval0.05_birVmal_BP.tsv', row.names = FALSE, quote=FALSE, sep="\t")
#write.table(go_output,'./TT-brain-33cV22c-w-mito_DGE_lfc-shr_pval0.05_33cV22cF1_BP.tsv', row.names = FALSE, quote=FALSE, sep="\t")
#write.table(go_output,'./TT-brain-33cV22c-w-mito_DGE_lfc-shr_pval0.05_33cV22cBir_BP.tsv', row.names = FALSE, quote=FALSE, sep="\t")
write.table(go_output,'./TT-brain-33cV22c-w-mito_DGE_lfc-shr_pval0.05_33cV22cMal_BP.tsv', row.names = FALSE, quote=FALSE, sep="\t")

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

#write.table(go_output,'./TT-brain-22c-w-mito_DGE_lfc-shr_pval0.05_birVmal_CC.tsv', row.names = FALSE, quote=FALSE, sep="\t")
#write.table(go_output,'./TT-brain-33c-w-mito_DGE_lfc-shr_pval0.05_birVmal_CC.tsv', row.names = FALSE, quote=FALSE, sep="\t")
#write.table(go_output,'./TT-brain-33cV22c-w-mito_DGE_lfc-shr_pval0.05_33cV22cF1_CC.tsv', row.names = FALSE, quote=FALSE, sep="\t")
#write.table(go_output,'./TT-brain-33cV22c-w-mito_DGE_lfc-shr_pval0.05_33cV22cBir_CC.tsv', row.names = FALSE, quote=FALSE, sep="\t")
write.table(go_output,'./TT-brain-33cV22c-w-mito_DGE_lfc-shr_pval0.05_33cV22cMal_CC.tsv', row.names = FALSE, quote=FALSE, sep="\t")


## random code below

# match sig go ids to gene nums
#go2gene_matchup <- data.frame("GOBPID" = genes_match$frame.go_id, "GeneID" = genes_match$frame.ensembl)

#c33high_go2gene_matchup <- genes_match
#write.csv(c33high_go2gene_matchup,'Data/dge/brain_w-mito/transgressive-exp/TT-brain-33c-w-mito_DGE_lfc-shr_all-F1trans.high-GOsiggenesmatch.csv')
# c33low_go2gene_matchup <- genes_match
# write.csv(c33low_go2gene_matchup,'Data/dge/brain_w-mito/transgressive-exp/TT-brain-33c-w-mito_DGE_lfc-shr_all-F1trans.low-GOsiggenesmatch.csv')

#c22high_go2gene_matchup <- genes_match
#write.csv(c22high_go2gene_matchup,'Data/dge/brain_w-mito/transgressive-exp/TT-brain-22c-w-mito_DGE_lfc-shr_all-F1trans.high-GOsiggenesmatch.csv')
c22low_go2gene_matchup <- genes_match
write.csv(c22low_go2gene_matchup,'Data/dge/brain_w-mito/transgressive-exp/TT-brain-22c-w-mito_DGE_lfc-shr_all-F1trans.low-GOsiggenesmatch.csv')

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
