## Approximate Bayesian Computation approach to estimate the
## range of effect sizes that may be explained by the chr22 QTL
## This script extracts observed genotypes and covariates at the QTL peak
## and assigns an expected CTmax based on genotype

## See Supporting Information 1 in Payne et al 2022 Supplementary Material

# load data file containing genotypes at all tested loci under the chromosome 22 QTL
data<-read.csv(file="Scripts/input_files/ltreb-ctmax-qtl_genos-chr22-qtl-peak.csv")
data<-na.omit(data)
stat1<-{}
stat2<-{}
r2<-{}
sim_effect<-{}

## 110k total simulations produces >500 accepted simulations
for(x in 1:110000){
  # randomly select a QTL effect size to simulate
  # ranging from 0% to 30%
  effect_size<-runif(1,0,0.3)
  Ct_diff<-3.4
  average_Ct<-36
  var_Ct<-0.42

  sim_effect<-c(sim_effect,effect_size)
  sim_data<-{}
  for(x in 1:length(data[,1])){

    if(data[,20][x] == 'AA' |data[,20][x] == 'BB' ){
      # if a parental/homozygous genotype under the chr22 QTL, a CTmax value is
      # randomly drawn from a normal distribution centered at the
      # average CTmax with standard deviation equal to the standard
      # deviation observed for the real CTmax values
      pheno<-rnorm(1,average_Ct,sqrt(var_Ct))
    }
    if(data[,20][x] == 'AB'){
      # if a hybrid/heterozygous genotype under the chr22 QTL, a CTmax value is
      # randomly drawn from a normal distribution centered at the
      # average CTmax with standard deviation equal to the standard
      # deviation observed for the real CTmax values, and then is subtracted
      # by the product of the simulated effect size of the QTL multiplied by
      # the average difference in CTmax between parental species
      pheno<-rnorm(1,average_Ct,sqrt(var_Ct))-effect_size*Ct_diff
    }

    sim_data<-c(sim_data,pheno)

  }

  ## run analysis
  # grab all relevant variables
  data_prob.df<-data          # dataset with simulated CTmax values
  ctmax <- sim_data           # simulated CTmax values
  hi <- data_prob.df$hi       # hybrid index
  # relevant covariates
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

  # grab the genotypes of the locus at the QTL peak
  geno <- data_prob.df$ScyDAA6.2113.HRSCAF.2539.8811359
  # organize variables of interest into new data frame, remove missing data
  dp_noNA<- data.frame(ctmax,hi,geno,STH.1,STH.2,STH.3,STH.4,STH.7,STH.8,STL.1,STL.2,STL.4,STL.5,STL.6,STL.7,STM.2,STM.4,STM.5,STM.6,STM.7)
  dp_noNA[dp_noNA == '-']<-NA
  dp_noNA<-na.omit(dp_noNA)
  #dim(dp_noNA)

  ## perform a linear regression with the full model and get the R^2 value
  lm.res<-lm(ctmax ~ geno + hi + STH.1 + STH.2 + STH.3 + STH.4 + STH.7 + STH.8 + STL.1 + STL.2 + STL.4 + STL.5 + STL.6 + STL.7 + STM.2 + STM.4 + STM.5 + STM.6 + STM.7,data=dp_noNA)
  full<-summary(lm.res)$adj.r.squared
  # perform a linear regression with the null model (no genotype) and get the R^2 value
  null.lm.res<-lm(ctmax ~ hi + STH.1 + STH.2 + STH.3 + STH.4 + STH.7 + STH.8 + STL.1 + STL.2 + STL.4 + STL.5 + STL.6 + STL.7 + STM.2 + STM.4 + STM.5 + STM.6 + STM.7,data=dp_noNA)
  null<-summary(null.lm.res)$adj.r.squared

 # record the difference in R^2 value between the full and null models
 # this is the estimated effect size for this simulation
  r2<-c(r2,full-null)
  # record additional statistics
  stat1<-c(stat1,summary(lm.res)$coefficients[,3][2]) # t-value
  stat2<-c(stat2,summary(lm.res)$coefficients[,4][2]) # Pr(>|t|) (p-value)
}

## accept or reject the simulation based on the estimated effect size, with a 5% tolerance threshold
## the estimated effect size r2 was 4.42%
threshold <- 0.05
observed_r2 <- 0.0442
r2_low  <- (1-threshold) * observed_r2
r2_high <- (1+threshold) * observed_r2

sim_stats <- as.data.frame(cbind(stat1,sim_effect,r2))
accepted<- subset(sim_stats,r2 <= r2_high & r2 >= r2_low)

## obtain 95% credibility interval of effect size consistent with our data
quantile(accepted$r2,c(0.025,0.975))
