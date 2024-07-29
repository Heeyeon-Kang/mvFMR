##############################################################################
####              Penalized estimation for a finite mixture of            ####
####                   multivariate regression models                     ####
##############################################################################
####                        Author: Heeyeon Kang                          ####
####                      Supervisor: Sunyoung Shin                       ####
##############################################################################
####       R-code of analyzing the result of the bike sharing data        ####
####                      presented in Section 6.2.                       ####
##############################################################################

# The following R-code analyzes the results of the bike sharing data and 
# generates Table 8.

# The results of the analysis of the bike sharing data are stored in 
# bike_sharing.rda, which consists of the list of output, density, w, and BIC.

# List of 4
# - output  : The modified BIC minimized by searching 39 lambdas 
#             for K = 1, 2, 3, 4, 5, 6.
# - w       : The conditional expectation of z_{ik} at E step which is 
#             the mixing proportion of the data points and is utilized 
#             for clustering.
# - density : The multivariate Gaussian distribution corresponding to each 
#             data point.
# - BIC     : The modified BIC for K = 1, 2, 3, 4, 5, 6.

load(file="./output/bike_sharing.rda")
source("./data/bike_sharing_data.R")
source("./code/functions.R")

optimal_K <- which.min(bike_sharing$BIC)
optimal_output <- bike_sharing$output[[optimal_K]]
optimal_w <- bike_sharing$w[[optimal_K]]
optimal_density <- bike_sharing$density[[optimal_K]]


### Generating the Table 7 ###
# Round down to five decimal places #
Bk_floor5 <- lapply(optimal_output$Bk, floor_5)

# To restore centering and scaling #
sd_X <- apply(X_before, 2, sd)
Bk_new <- Bk_floor5
for(k in 1:optimal_K){
  for(i in 10:ncol(X)){
    for(j in 1:ncol(Y)){
      if(Bk_floor5[[k]][i,j] != 0){
        Bk_new[[k]][i,j] <- Bk_floor5[[k]][i,j] / sd_X[i]
      }
    }
  }
}

# Calculating the correlation between glucose and HbA1c #
correlation <- sapply(optimal_output$sigma, function(x) x[1,2]/sqrt(x[1,1]*x[2,2]))

# Combining to make the Table 8 #
bike_sharing_analysis <- list(pi=optimal_output$pi, Bk=Bk_new, corr=correlation)


