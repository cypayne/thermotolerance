##  Perform NULL admixture mapping with a genotypes file, a hybrid index file, 
##  a phenotypes file, a column number for the focal phenotype of interest, 
##  a start column index marking where covariate columns start, an end column
##  index marking where covariate columns end, and an outfile name tag to 
##  be appended to the outfile name
##  NOTE: sample order *must* be identical across files
##  This code is similar to that in perform_glm_admixture_mapping_v2_gaussian.R, 
##  but with additional covariates

### WARNING: This script was written for a specific case and has hardcoded elements, 
###          needs to be edited if not using 17 site-tank covariates

### NULL simulations for perform_glm_admixture_mapping_v2_gaussian_plus-multiple-cov.R
### usage: Rscript perform_glm_admixture_mapping_v2_gaussian_plus-multiple-cov.R 
###                genotypes_file hybrid_index_file phenotypes_file focal_column_number 
###                covar_column_start_number covar_column_end_number name_tag

### ms & cyp 06/2020

# To get helpful output for debugging, set to TRUE
DEBUG = FALSE 

# Check that correct number of arguments are passed
args <- commandArgs(TRUE)
if(length(args)<7){
stop("usage is Rscript perform_glm_admixture_mapping_v2_gaussian_plus-one-cov.R genotypes_file phenotypes_file focal_column_number covar_column_start_number covar_column_end_number name_tag num_sims")
}

# Load arguments
infile        <- as.character(args[1])
data          <-read.csv(file=infile,head=T,as.is=T,sep="\t")
pheno         <-as.character(args[2])
phenotypes    <-read.csv(file=pheno,sep="\t",head=TRUE)
pheno_column  <-as.numeric(args[3])
covar_start_column  <-as.numeric(args[4])
covar_end_column    <-as.numeric(args[5])
tag           <-as.character(args[6])
num_sims      <-as.numeric(args[7])

# Run through each simulation
for(y in 1:num_sims){
  out<-paste(infile,"_results_gaussian_",tag,"_",y,sep="")
  file.remove(out)

  names<-colnames(data)
  combined_phenos<-cbind(phenotypes[,pheno_column],phenotypes[,covar_start_column:covar_end_column])
  combined_phenos_shuff<-combined_phenos[sample(nrow(combined_phenos),length(combined_phenos[,1])),] 
  null_phenos<-combined_phenos_shuff[,1]
  null_covars<-combined_phenos_shuff[,2:length(combined_phenos_shuff)]

  # helpful for debugging
  if(y==1 | y==2) {
    print(head(combined_phenos))
    print(head(combined_phenos_shuff))
    print(head(null_covars))
  }

  ## Perform admixture mapping for each marker
  track=0
  for (x in 2:length(data[1,])){
    cat("sim ", y, ", marker ", x,"\n")

    # collect pheno column, covar column, hi column, and marker genotypes
    dat<-cbind(null_phenos,phenotypes$hybrid_index,null_covars,data[,x])
    last_index<-length(dat)

    # omit individuals with NAs 
    dat<-na.omit(dat)

    # remove marker if less than 1/2 of the individuals are genotyped
    if(length(dat[,1]) >= (length(data[,1])*0.5) & length(unique(dat[,last_index]))>1){

      # run null model (i.e. pheno~index+covariates)
      # having an error issue here - contrasts error pops up for some
      # markers, not sure why - need to throw those out
      potentialError<- tryCatch(
        model2<-glm(as.numeric(dat[,1])~as.numeric(dat[,2]) + as.factor(dat[,3]) + 
                  as.factor(dat[,4]) + as.factor(dat[,5]) +
                  as.factor(dat[,6]) + as.factor(dat[,7]) + as.factor(dat[,8]) + 
                  as.factor(dat[,9]) + as.factor(dat[,10]) +
                  as.factor(dat[,11]) + as.factor(dat[,12]) + as.factor(dat[,13]) + 
                  as.factor(dat[,14]) + as.factor(dat[,15]) + as.factor(dat[,16]) + 
                  as.factor(dat[,17]) + as.factor(dat[,18]) +
                  as.factor(dat[,19]),family="gaussian"), error=function(e) e )

      if(inherits(potentialError, "error")) {
        cat("threw contrasts error\n")
        next
      }

      null<-logLik(model2)[1]

      # run full model (i.e. pheno~index+covariates+genotype
      # genotype coded as factor here
      model1<-glm(as.numeric(dat[,1])~as.numeric(dat[,2]) + as.factor(dat[,3]) + 
                  as.factor(dat[,4]) + as.factor(dat[,5]) +
                  as.factor(dat[,6]) + as.factor(dat[,7]) + as.factor(dat[,8]) + 
                  as.factor(dat[,9]) + as.factor(dat[,10]) +
                  as.factor(dat[,11]) + as.factor(dat[,12]) + as.factor(dat[,13]) + 
                  as.factor(dat[,14]) + as.factor(dat[,15]) + as.factor(dat[,16]) + 
                  as.factor(dat[,17]) + as.factor(dat[,18]) +
                  as.factor(dat[,19]) + as.factor(dat[,20]),family="gaussian")

      focal<-logLik(model1)[1]

      # get log likelihood difference
      like_diff<-focal-null

      # grab intercept, hi, and marker p-vals
      p_list <- summary(model1)$coef[, "Pr(>|t|)"]
      if(y==1 | y==2) { print(summary(model1)) }
      results<-cbind(names[x],p_list[1],p_list[2],p_list[length(p_list)-1],p_list[length(p_list)],(summary(model1)$coef[,"t value"])[length(p_list)],like_diff,length(dat[,1]))

      # cat((summary(model1)$coef[,"t value"])[length(p_list)],"\n",like_diff,"\n")

      # Write out results (one p-value per categorical genotype) 
      if(track==0){
        write.table(results,file=out,append=TRUE,col.names=c("chrom.marker","intercept","mixture_prop","geno1_pval","geno2_pval","t-value","likelihood-diff","num_ind"),row.names=F,sep="\t",quote=FALSE)
        track=1
      } else{
        write.table(cbind(names[x],p_list[1],p_list[2],p_list[length(p_list)-1],p_list[length(p_list)],(summary(model1)$coef[,"t value"])[length(p_list)],like_diff,length(dat[,1])),file=out,append=TRUE,col.names=FALSE,row.names=F,sep="\t",quote=FALSE)
      }

    }#if half the data

  }#for all lines

} #for all sims
