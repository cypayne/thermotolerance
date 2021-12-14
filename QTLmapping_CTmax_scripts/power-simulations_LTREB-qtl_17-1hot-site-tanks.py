#!/usr/bin/env python3

''' Determine likelihood of effect size for a qtl
    using power simulations (power threshold)
    i.e. how frequently do you recover a qtl of 
    this estimated effect size by chance, how much
    should we trust qtl effect

    Simulates phenotypes for observed genotypes under qtl

    This script is for LTREB ctmax qtl

    individual_file: contains column "hybrid_index" per
                     observed individual, as well as 
                     relevant covariates
    genotype_file:   contains single column called 
                     "genotype" with observed genos
                     per indv under qtl peak
    cov_start_index: [if including covariates]
                     first covariate index
    cov_end_index:   [if including covariates]
                     last covariate index
    
    usage: ./power-simulations.py individual_file genotype_file [cov_start_index cov_end_index]

    for a underdominant qtl:
      heterozygotes should have qtl_effect lower commpared to parents 

    cyp xi-2020
'''

import sys
import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf

if len(sys.argv) < 3:
  print("Provide a file with 'hybrid_index' column")
  exit()
elif len(sys.argv) == 5: 
  cov_start_index = int(sys.argv[3])
  cov_end_index = int(sys.argv[4])

num_sims     = 10 # number of simulations

qtl_effect   = 0.0442    # effect size, variation explained
#pheno_total  = 39.2   # avg ctmax of birchmanni (BB)
pheno_total  = 36.0   # avg ctmax of indvs with homo genotype (BB and MM) 
pheno_qtl    = (pheno_total-pheno_other)*qtl_effect
pheno_other  =        # avg ctmax of indvs with het genotype (BM)
stdev        = 1.026517 # trait stdev
## use corrected p-value, threshold for gw significance
## computed from null simulations (most min p-val)
pval_cutoff  = 0.0000244 # 95% p-value accepted
pval_list    = []

# dictionary of genotype to phenotype 
#geno2pheno_dict = \
#  { 0: pheno_other,
#    1: pheno_other + pheno_qtl*0.5,
#    2: pheno_total
#  }

# underdominant qtl
geno2pheno_dict = \
  { 0: pheno_total,
    1: pheno_other,
    2: pheno_total
  } 

# pass in file with all observed indvs, and at least
# one col with hybrid_index per individual
indv_df = pd.read_csv(sys.argv[1],sep="\t")
indv_ancestries = indv_df["hybrid_index"].tolist()
#for cov_index in range(cov_start_index,cov_end_index):
#  covariate_list = indv_df.iloc[:,cov_index].tolist()

# set seed
np.random.seed(84123)

for sim in range(0,num_sims):
  print("sim {}...".format(sim))
  geno_vector, pheno_vector = [], []
  # randomly simulate genotypes, based on average 
  # ancestry for many individuals

  # ancestry for each indv as observed + tanks as observed
  for avg_ancestry in indv_ancestries:
    # generate genotype per chrom, around avg ancestry
    try:
      rand_geno = np.random.binomial(1,avg_ancestry,2)
    except Exception as e:
      print(avg_ancestry, type(avg_ancestry))
      raise e
    # get genotype code by adding two genotypes together
    geno = rand_geno[0]+rand_geno[1]
    geno_vector.append(geno)

    # add trait variation not explained by qtl
    # mean is equal to pheno_other, stdev = observed trait stdev in F2s
    rand_pheno_var = np.random.normal(pheno_other, stdev, 1)[0] 
    pheno = geno2pheno_dict[geno] + (rand_pheno_var - pheno_other)
    pheno_vector.append(pheno)

  # and then look for qtl
  # create dataframe from lists
  dat = pd.DataFrame({'phenotype': pheno_vector, 'genotype': geno_vector, \
                     'hybrid_index': indv_ancestries})

  formula = 'phenotype ~ hybrid_index + genotype'
  #null_formula = 'phenotype ~ hybrid_index'

  if cov_start_index: 
    # add covariate columns to dat df
    for cov_index in range(cov_start_index,cov_end_index+1):
      # wrap the name in C() to denote cov entry as a categorical var
      cov_name = 'cov'+str(cov_index)
      dat[cov_name] = indv_df.iloc[:,cov_index] 

    # create formula based on covariates
    #formula = 'phenotype ~ hybrid_index + C(covariate) + genotype'
    #null_formula = 'phenotype ~ hybrid_index + C(covariate)'
    formula = 'phenotype ~ C(genotype) + hybrid_index + ' + \
              ' + '.join(['C('+col+')' for col in dat.columns if not col=="phenotype" and not col=="genotype" and not col=="hybrid_index"])

  print(dat)
  print(formula)
  #null  = smf.glm(formula=null_formula,data=dat,family=sm.families.Gaussian())
  model = smf.glm(formula=formula,data=dat,family=sm.families.Gaussian())
  result = model.fit()

  print(result.pvalues)
  ## not sure if this intuition is right, but pull the pvalue for the heterozygous genotype (since this is typically comparable to pvalue
  ## when genotype is numerically encoded
  ## i.e. collect pval after intercept pval, so position 1 --> C(genotype)[T.1]
  pval = result.pvalues[1]

  pval_list.append(pval)

num_times_qtl_recovered = sum(i < pval_cutoff for i in pval_list)
print(pval_list)
print("Parameter summary: num_sims={}, qtl_effect={}, pheno_total={}, stdev={}, cutoff={}, formula={}, infile={}".format(num_sims,qtl_effect,pheno_total,stdev,pval_cutoff, formula,sys.argv[1]))
print("Of {} simulations, {} qtl were recovered.".format(num_sims, num_times_qtl_recovered))
