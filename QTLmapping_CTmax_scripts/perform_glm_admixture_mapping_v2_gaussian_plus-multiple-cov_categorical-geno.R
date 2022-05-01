## perform_glm_admixture_mapping_v2_gaussian_plus-multiple-cov_categorical-geno.R
#
##  Perform admixture mapping with a genotypes file, a hybrid index file,
##  a phenotypes file, a column number for the focal phenotype of interest,
##  a start column index marking where covariate columns start, an end column
##  index marking where covariate columns end, and an outfile name tag to
##  be appended to the outfile name
##  NOTE: sample order *must* be identical across files
##  This code is similar to that in perform_glm_admixture_mapping_v2_gaussian.R,
##  but with additional covariates

### WARNING: This script was written for a specific case and has hardcoded elements,
###          needs to be edited if not using 17 site-tank covariates

### Essentially the same as perform_glm_admixture_mapping_v2_gaussian_plus-multiple-cov.R,
### but with genotype encoded as a categorical variable (i.e. 3 categories)
### usage: Rscript perform_glm_admixture_mapping_v2_gaussian_plus-one-cov.R
###                genotypes_file hybrid_index_file phenotypes_file focal_column_number
###                covar_column_start_number covar_column_end_number name_tag

### ms & cyp 06/2020

# To get helpful output for debugging, set to TRUE
DEBUG = FALSE

# Check that correct number of arguments are passed
args <- commandArgs(TRUE)
if(length(args)<7){
  stop("usage is Rscript perform_glm_admixture_mapping_v2_gaussian_plus-one-cov.R genotypes_file hybrid_index_file phenotypes_file focal_column_number covar_column_start_number covar_column_end_number name_tag")
}

# Load arguments
infile        <- as.character(args[1])
data          <-read.csv(file=infile,head=T,as.is=T,sep="\t")
hybrid_index  <-as.character(args[2])
index         <-read.csv(file=hybrid_index,sep="\t",head=TRUE)
pheno         <-as.character(args[3])
phenotypes    <-read.csv(file=pheno,sep="\t",head=TRUE)
pheno_column  <-as.numeric(args[4])
covar_start_column  <-as.numeric(args[5])
covar_end_column    <-as.numeric(args[6])
tag           <-as.character(args[7])
out           <-paste(infile,"_results_gaussian_v2_",tag,sep="")
file.remove(out)

names<-colnames(data)

## Perform admixture mapping for each marker
track=0
for (x in 2:length(data[1,])){
  # collect pheno column, covar column, hi column, and marker genotypes
  dat<-cbind(phenotypes[,pheno_column],index$hybrid_index,phenotypes[,covar_start_column:covar_end_column],data[,x])
  last_index<-length(dat)

  # omit individuals with NAs
  dat<-na.omit(dat)

  # helpful for debugging
  if(DEBUG) {
    if(x==2) {
      print(head(dat))
      print(dim(dat))
    }
  }

  # remove marker if less than 1/2 of the individuals are genotyped
  if(length(dat[,1]) >= (length(data[,1])*0.5) & length(unique(dat[,last_index]))>1){

    # run null model (i.e. pheno~index)
    model2<-glm(as.numeric(dat[,1])~as.numeric(dat[,2]) + as.factor(dat[,3]) +
                as.factor(dat[,4]) + as.factor(dat[,5]) +
                as.factor(dat[,6]) + as.factor(dat[,7]) + as.factor(dat[,8]) +
                as.factor(dat[,9]) + as.factor(dat[,10]) +
                as.factor(dat[,11]) + as.factor(dat[,12]) + as.factor(dat[,13]) +
                as.factor(dat[,14]) + as.factor(dat[,15]) + as.factor(dat[,16]) +
                as.factor(dat[,17]) + as.factor(dat[,18]) +
                as.factor(dat[,19]), family="gaussian")

    null<-logLik(model2)[1]

    # run full model (i.e. pheno~index+covariates+genotype
    # genotype coded as numeric here
    model1<-glm(as.numeric(dat[,1])~as.numeric(dat[,2]) + as.factor(dat[,3]) +
                as.factor(dat[,4]) + as.factor(dat[,5]) +
                as.factor(dat[,6]) + as.factor(dat[,7]) + as.factor(dat[,8]) +
                as.factor(dat[,9]) + as.factor(dat[,10]) +
                as.factor(dat[,11]) + as.factor(dat[,12]) + as.factor(dat[,13]) +
                as.factor(dat[,14]) + as.factor(dat[,15]) +as.factor(dat[,16]) +
                as.factor(dat[,17]) + as.factor(dat[,18]) + as.factor(dat[,19]) +
                as.factor(dat[,20]),family="gaussian")

    focal<-logLik(model1)[1]

    # get log likelihood difference
    like_diff<-focal-null

    # grab intercept, hi, and marker p-vals
    p_list <- summary(model1)$coef[, "Pr(>|t|)"]
    print(summary(model1))
    results<-cbind(names[x],p_list[1],p_list[2],p_list[length(p_list)-1],p_list[length(p_list)],(summary(model1)$coef[,"t value"])[length(p_list)],like_diff,length(dat[,1]))

    if(DEBUG) {
      cat((summary(model1)$coef[,"t value"])[length(p_list)],"\n",like_diff,"\n")
      cat(results)
    }

    # Write out results
    if(track==0){
      write.table(results,file=out,append=TRUE,col.names=c("chrom.marker","intercept","mixture_prop","geno1_pval","geno2_pval","t-value","likelihood-diff","num_ind"),row.names=F,sep="\t",quote=FALSE)
      track=1
    } else{
      write.table(cbind(names[x],p_list[1],p_list[2],p_list[length(p_list)-1],p_list[length(p_list)],(summary(model1)$coef[,"t value"])[length(p_list)],like_diff,length(dat[,1])),file=out,append=TRUE,col.names=FALSE,row.names=F,sep="\t",quote=FALSE)

    }

  }#if half the data

}#for all lines
