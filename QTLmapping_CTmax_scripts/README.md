# Pipeline for mapping CTmax QTL 

This is a cleaned pipeline for mapping QTL with both R/qtl and a general linear model (GLM). 
This analysis is presented in Payne et al 2021, where we mapped the critical thermal 
maximum (CTmax) of ~150 Xiphophorus malinche-X. birchmanni F2 hybrids. Below is a
summary of the included scripts.

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

                   
* **power-simulations_LTREB-qtl_17-1hot-site-tanks.py** 

  * *Description: Perform power simulations to determine the likelihood of recovering a 
                  QTL of some effect size by chance.*


## Contact

Questions and suggested improvements are greatly appreciated - please reach out to 
the repository owner at cypayne at stanford dot edu.
