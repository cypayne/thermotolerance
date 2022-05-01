# QTL mapping 

This is a clean pipeline for mapping QTL with both R/qtl and a general linear model (GLM). 
This analysis is presented in Payne et al 2022, where we mapped the critical thermal 
maximum (CTmax) of ~150 Xiphophorus malinche-X. birchmanni artificial F2 hybrids.

All input files used in the following analysis can be found in input_files/ .

## Dependencies

R-3.6.1

The following packages are required for each anlaysis (see install.R):
* QTL-mapping ~ qtl

## Contents

* **CTmax-Rqtl-analysis_clean.R**

  * *Description: Full R/qtl pipeline, including sample and marker filtration steps,
                  model selection, QTL scan, LOD threshold permutation, and plots.*

* **perform_glm_admixture_mapping_v2_gaussian_plus-multiple-cov.R** 

  * *Description: Perform admixture mapping (with GLM), where genotype is encoded
                  as a numeric variable.* 

* **perform_glm_admixture_mapping_v2_gaussian_plus-multiple-cov_categorical-geno.R** 

  * *Description: Perform admixture mapping (with GLM), where genotype is encoded
                  as a categorical variable.*

* **perform_glm_admixture_mapping_v2_gaussian_plus-multiple-cov_categorical-geno_NULL.R**
  
  * *Description: Perform null simulations of admixture mapping (categorical genotype) to
                  get permuted 5% p-value GLM threshold.*

* **CTmax-GLM-analysis_clean.R**

  * *Description: Plot results of admixture mapping.* 

* **ABC_simulations_for_effect_size.R**

  * *Description: Approximate Bayesian Computation approach to estimate the
                  range of effect sizes that may be explained by the QTL.*

* **power-simulations_LTREB-qtl_17-1hot-site-tanks.py** 

  * *Description: Perform power simulations to determine the likelihood of recovering a 
                  QTL of some effect size by chance.*


## Contact

Questions and suggested improvements are greatly appreciated - please reach out to 
the repository owner at cypayne at stanford dot edu.
