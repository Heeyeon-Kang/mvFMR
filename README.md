# mvFMR

### Electronic supplement to "Penalized estimation for a finite mixture of regression models"

This repository has the R-code to do the analyses and reproduce the figures and tables in Sections 5 and Section 6 of the manuscript and the supplementary material.

To run the R-code, it is recommended to load the R project "FMRwithMultipleResponses.Rproj", as all paths are set relative to this directory.

The repository consists of the following folders:
* Data: R-code for generating or refining the data used in Section 5 and Section 6.
  * "simulation_data.R" contains the functions generating the dataset using in Section 5. 
* Functions: R-code of all functions for running the EM-ADMM algorithm and R-code for fitting simulation data using each method.
* Outputs: The results of simulation studies and real data analyses to reproduce the figures and tables presented in Section 5 and Section 6
