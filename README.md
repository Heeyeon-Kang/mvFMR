# mvFMR

### Electronic supplement to "Penalized estimation for a finite mixture of regression models"


This repository has the R-code to do the analyses and reproduce the figures and tables in Sections 5 and Section 6 of the manuscript and the supplementary material.

All code was written in a laptop with R version 4.3.1 and run on a Linux Cluster with R version 4.1.2.
To run the R-code, it is recommended to load the R project "FMRwithMultipleResponses.Rproj", as all paths are set relative to this directory.

The repository consists of the following folders:

* Data: R-code for generating or refining the data used in Section 5 and Section 6;
  * "simulation_seed_number.R" contains the seed numbers generating the data of simulations in Section 5.
  * <span style="color:orange"> simulation_data.R </span> contains the functions.
  * "simulation_data.R" contains the functions generating the dataset using in Section 5.
  * "diabetes_diagnosis_data.R" and "bike_sharing_data.R" are the R-codes of the process of refining the raw data, "diabetes_diagnosis.csv" and "bike_sharing.csv", respectively.
 
* Functions: R-code of all functions for running the EM-ADMM algorithm and R-code for fitting simulation data using each method;
  * "functions.R" is the R-code of all functions for running EM-ADMM algorithm.
  * “mvFMR.R”, "mvFMR_LASSO.R", "mvFMR_SCAD.R", and "mvFMR_MCP.R" contain the R-codes for fitting simulation data using mvFMR with their respective penalty functions presented in Section 5.3.
  * "mvFMR_fixedK.R", "mvFMR_LASSO_fixedK.R", "mvFMR_SCAD_fixedK.R", and "mvFMR_MCP_fixedK.R" contain the R-codes for fitting the data using mvFMR with their respective penalty functions presented in Section 5.1, Section 5.2, and Section 5.3, assuming K is known.
  * "flexmix_fixedK.R" contains the R-code for fitting the data using the R-package "Flexmix" presented in Section 5.1 and Section 5.2.
  * "oracle_fixedK.R" contains the R-code to estimate the oracle estimator presented in Section 5.1 and Section 5.2 by fitting the data.

* Simulations: The results of simulation studies and real data analyses to reproduce the figures and tables presented in Section 5 and Section 6;
  * "simulation_table.R" contains the functions for calculating TPR, FPR, MSE, and predictive log-likelihood loss.
  * "diabetes_diagnosis_analysis_mvFMR_S.R" and "bike_sharing_analysis_mvFMR_S.R" contain the R-codes for analyzing the results of the data and generating Table 7 and Table 8 presented in Section 6.1 and Section 6.2.
  * "diabetes_diagnosis_analysis_k-means.R" contains the R-codes for analyzing the result of the diabetes diagnosis data by using k-means clustering and generating Figure 1 presented in Section 6.1.
  * "diabetes_diagnosis_fitting.R" and "bike_sharing_fitting.R" contain the R-codes for fitting the data presented in Section 6.1 and Section 6.2 using parallel computing.
  * "s1_table.rda", "s2_table.rda", and "s3_table.rda" are the rda files containing the contents of the tables in Section 5.1, Section 5.2, and Section 5.3, respectively.
  * "diabetes_diagnosis.rda" and "bike_sharing.rda" are rda files containing the results of the analysis using mvFMR-SCAD.

    
