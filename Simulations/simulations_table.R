##############################################################################
####              Penalized estimation for a finite mixture of            ####
####                   multivariate regression models                     ####
##############################################################################
####                        Author: Heeyeon Kang                          ####
####                      Supervisor: Sunyoung Shin                       ####
##############################################################################
####          R-code for the descriptions for Table 2, 4, and 6.          ####
##############################################################################

# The following R-code contains functions for calculating TPR, FPR, MSE, and 
# predictive log-likelihood loss.

# "s1_table.rda" is the rda file of the outputs in Section 5.1.
# The s1_table consists of model1, model2, and model3, with the mean and sd 
# for each model if the sample size is 500, 1000, or 2000.
# "mean" and "sd" are the mean and standard deviation of TPR, FPR, and MSE of
# 500 independent simulated trials.
# List of 3
# - model1 
#    - "500" : size n=500 in model 1
#       - mean
#       - sd
#    - "1000" : size n=1000 in model 1
#       - mean
#       - sd
# - model2
#    - "500" : size n=500 in model 2
#       - mean
#       - sd
#    - "1000" : size n=1000 in model 2
#       - mean
#       - sd
# - model3
#    - "1000" : size n=1000 in model 3
#       - mean
#       - sd
#    - "2000" : size n=2000 in model 3
#       - mean
#       - sd

# "s2_table.rda" is the rda file of the outputs in Section 5.2.
# The s2_table consists of model1, model2, and model3, with the mean and sd
# for each model if the sample size is 500 and 1000.
# "mean" and "sd" are the mean and standard deviation of TPR, FPR, and MSE of
# 500 independent simulated trials.
# List of 3
# - model1 
#    - "500" : size n=500 in model 1
#       - mean
#       - sd
#    - "1000" : size n=1000 in model 1
#       - mean
#       - sd
# - model2
#    - "500" : size n=500 in model 2
#       - mean
#       - sd
#    - "1000" : size n=1000 in model 2
#       - mean
#       - sd
# - model3
#    - "500" : size n=500 in model 3
#       - mean
#       - sd
#    - "1000" : size n=1000 in model 3
#       - mean
#       - sd

# "s3_table.rda" is the rda file of the outputs in Section 5.3.
# The s3_table consists of model1 and model2 with the mean and sd for each model
# if the sample size is 500 and 1000.
# "mean" and "sd" are the mean and standard deviation of predictive log-likelihood loss
# of 500 independent simulated trials.
# List of 2
# - model1 
#    - "500" : size n=500 in model 1
#       - mean
#       - sd
#    - "1000" : size n=1000 in model 1
#       - mean
#       - sd
# - model2
#    - "500" : size n=500 in model 2
#       - mean
#       - sd
#    - "1000" : size n=1000 in model 2
#       - mean
#       - sd

TPR <- function(true, pred, decimal=c("5","6")){
  
  total_true_Bk <- sapply(true$Bk, as.vector)
  total_pred_Bk <- sapply(pred$Bk, as.vector)
  
  Number_of_nonzeros <- length(which(total_true_Bk != 0))
  
  norm_of_true_and_pred_Bk <- vector()
  for(i in 1:nrow(permutations(K,K))){
    norm_of_true_and_pred_Bk[i] <- norm(total_true_Bk - total_pred_Bk[, permutations(K,K)[i,]], "1")
  }
  switch_of_pred <- permutations(K,K)[which.min(norm_of_true_and_pred_Bk),]
  
  for(j in 2:4){
    pred[[j]] <- pred[[j]][switch_of_pred]
  }
  
  total_pred_Bk <- sapply(pred$Bk, as.vector)
  
  if(decimal == "5"){
    True_Positive_Rate <- length(which(total_true_Bk != 0 & floor_5(total_pred_Bk) != 0)) / Number_of_nonzeros
  }else if(decimal == "6"){
    True_Positive_Rate <- length(which(total_true_Bk != 0 & floor_6(total_pred_Bk) != 0)) / Number_of_nonzeros
  }else{
    message("Please select 5 or 6")
  }
  return(True_Positive_Rate)
}
FPR <- function(true, pred, decimal=c("5","6")){
  
  total_true_Bk <- sapply(true$Bk, as.vector)
  total_pred_Bk <- sapply(pred$Bk, as.vector)
  
  Number_of_zeros <- length(which(total_true_Bk == 0))
  
  norm_of_true_and_pred_Bk <- vector()
  for(i in 1:nrow(permutations(K,K))){
    norm_of_true_and_pred_Bk[i] <- norm(total_true_Bk - total_pred_Bk[, permutations(K,K)[i,]], "1")
  }
  switch_of_pred <- permutations(K,K)[which.min(norm_of_true_and_pred_Bk),]
  
  for(j in 2:4){
    pred[[j]] <- pred[[j]][switch_of_pred]
  }
  
  total_pred_Bk <- sapply(pred$Bk, as.vector)
  
  if(decimal == "5"){
    False_Positive_Rate <- length(which(total_true_Bk == 0 & floor_5(total_pred_Bk) != 0)) / Number_of_zeros
  }else if(decimal == "6"){
    False_Positive_Rate <- length(which(total_true_Bk == 0 & floor_6(total_pred_Bk) != 0)) / Number_of_zeros
  }else{
    message("Please select 5 or 6")
  }
  return(False_Positive_Rate)
}
TPR_no_penalty <- function(true, pred){
  
  total_true_Bk <- sapply(true$Bk, as.vector)
  total_pred_Bk <- sapply(pred$Bk, as.vector)
  
  Number_of_nonzeros <- length(which(total_true_Bk != 0))
  
  norm_of_true_and_pred_Bk <- vector()
  for(i in 1:nrow(permutations(K,K))){
    norm_of_true_and_pred_Bk[i] <- norm(total_true_Bk - total_pred_Bk[, permutations(K,K)[i,]], "1")
  }
  switch_of_pred <- permutations(K,K)[which.min(norm_of_true_and_pred_Bk),]
  
  for(j in 1:3){
    pred[[j]] <- pred[[j]][switch_of_pred]
  }
  
  total_pred_Bk <- sapply(pred$Bk, as.vector)
  
  True_Positive_Rate <- length(which(total_true_Bk != 0 & total_pred_Bk != 0)) / Number_of_nonzeros
  
  return(True_Positive_Rate)
}
FPR_no_penalty <- function(true, pred){
  
  
  total_true_Bk <- sapply(true$Bk, as.vector)
  total_pred_Bk <- sapply(pred$Bk, as.vector)
  
  Number_of_zeros <- length(which(total_true_Bk == 0))
  
  norm_of_true_and_pred_Bk <- vector()
  for(i in 1:nrow(permutations(K,K))){
    norm_of_true_and_pred_Bk[i] <- norm(total_true_Bk - total_pred_Bk[, permutations(K,K)[i,]], "1")
  }
  switch_of_pred <- permutations(K,K)[which.min(norm_of_true_and_pred_Bk),]
  
  for(j in 1:3){
    pred[[j]] <- pred[[j]][switch_of_pred]
  }
  
  total_pred_Bk <- sapply(pred$Bk, as.vector)
  
  False_Positive_Rate <- length(which(total_true_Bk == 0 & total_pred_Bk != 0)) / Number_of_zeros
  
  return(False_Positive_Rate)
}
TPR_flexmix <- function(true, pred, decimal=c("5","6")){
  
  total_true_Bk <- sapply(true$Bk, as.vector)
  total_pred_Bk <- sapply(pred$Bk, as.vector)
  
  Number_of_nonzeros <- length(which(total_true_Bk != 0))
  
  norm_of_true_and_pred_Bk <- vector()
  for(i in 1:nrow(permutations(K,K))){
    norm_of_true_and_pred_Bk[i] <- norm(total_true_Bk - total_pred_Bk[, permutations(K,K)[i,]], "1")
  }
  switch_of_pred <- permutations(K,K)[which.min(norm_of_true_and_pred_Bk),]
  
  #pred <- pred[switch_of_pred]
  for(j in 1:3){
    pred[[j]] <- pred[[j]][switch_of_pred]
  }
  
  total_pred_Bk <- sapply(pred$Bk, as.vector)
  
  if(decimal == "5"){
    True_Positive_Rate <- length(which(total_true_Bk != 0 & floor_5(total_pred_Bk) != 0)) / Number_of_nonzeros
  }else if(decimal == "6"){
    True_Positive_Rate <- length(which(total_true_Bk != 0 & floor_6(total_pred_Bk) != 0)) / Number_of_nonzeros
  }else{
    message("Please select 5 or 6")
  }
  return(True_Positive_Rate)
}
FPR_flexmix <- function(true, pred, decimal=c("5","6")){
  total_true_Bk <- sapply(true$Bk, as.vector)
  total_pred_Bk <- sapply(pred$Bk, as.vector)
  
  Number_of_zeros <- length(which(total_true_Bk == 0))
  
  norm_of_true_and_pred_Bk <- vector()
  for(i in 1:nrow(permutations(K,K))){
    norm_of_true_and_pred_Bk[i] <- norm(total_true_Bk - total_pred_Bk[, permutations(K,K)[i,]], "1")
  }
  switch_of_pred <- permutations(K,K)[which.min(norm_of_true_and_pred_Bk),]
  
  for(j in 1:3){
    pred[[j]] <- pred[[j]][switch_of_pred]
  }
  
  total_pred_Bk <- sapply(pred$Bk, as.vector)
  
  if(decimal == "5"){
    False_Positive_Rate <- length(which(total_true_Bk == 0 & floor_5(total_pred_Bk) != 0)) / Number_of_zeros
  }else if(decimal == "6"){
    False_Positive_Rate <- length(which(total_true_Bk == 0 & floor_6(total_pred_Bk) != 0)) / Number_of_zeros
  }else{
    message("Please select 5 or 6")
  }
  return(False_Positive_Rate)
}
MSE <- function(true, pred){
  
  total_true_Bk <- sapply(true$Bk, as.vector)
  total_pred_Bk <- sapply(pred$Bk, as.vector)
  
  norm_of_true_and_pred_Bk <- vector()
  for(i in 1:nrow(permutations(K,K))){
    norm_of_true_and_pred_Bk[i] <- norm(total_true_Bk - total_pred_Bk[, permutations(K,K)[i,]], "1")
  }
  switch_of_pred <- permutations(K,K)[which.min(norm_of_true_and_pred_Bk),]
  
  for(j in 2:4){
    pred[[j]] <- pred[[j]][switch_of_pred]
  }
  
  total_pred_Bk <- sapply(pred$Bk, as.vector)
  total_pred_sigma <- sapply(pred$sigma, as.vector)
  total_true_sigma <- sapply(true$sigma, as.vector)
  
  MSE_pi <- 1/K * (norm(true$pi-pred$pi, "2"))^2
  MSE_Bk <- 1/K * (norm(total_true_Bk-total_pred_Bk, "F"))^2
  MSE_sigma <- 1/K * (norm(total_true_sigma-total_pred_sigma, "F"))^2
  
  return(c(MSE_pi, MSE_Bk, MSE_sigma))
}
MSE_no_penalty <- function(true, pred){
  
  total_true_Bk <- sapply(true$Bk, as.vector)
  total_pred_Bk <- sapply(pred$Bk, as.vector)
  
  norm_of_true_and_pred_Bk <- vector()
  for(i in 1:nrow(permutations(K,K))){
    norm_of_true_and_pred_Bk[i] <- norm(total_true_Bk - total_pred_Bk[, permutations(K,K)[i,]], "1")
  }
  switch_of_pred <- permutations(K,K)[which.min(norm_of_true_and_pred_Bk),]
  
  for(j in 1:3){
    pred[[j]] <- pred[[j]][switch_of_pred]
  }
  
  total_pred_Bk <- sapply(pred$Bk, as.vector)
  total_pred_sigma <- sapply(pred$sigma, as.vector)
  total_true_sigma <- sapply(true$sigma, as.vector)
  
  MSE_pi <- 1/K * (norm(true$pi-pred$pi, "2"))^2
  MSE_Bk <- 1/K * (norm(total_true_Bk-total_pred_Bk, "F"))^2
  MSE_sigma <- 1/K *(norm(total_true_sigma-total_pred_sigma, "F"))^2
  
  return(c(MSE_pi, MSE_Bk, MSE_sigma))
}
pred_ll <- function(l, n, output, K, w, density){
  
  # n : the size of data
  optimal_K <- K[[l]]
  w_ <- w[[l]]
  pi_ <- output[[l]]$pi
  dens <- density[[l]]
  
  logpi_term <- sum(sapply(w_, sum) * log(pi_))
  logdensity_term <- 0
  for(k in 1:optimal_K){
    logdensity_term <- logdensity_term + sum(w_[[k]] * lapply(dens, log)[[k]])
  }
  pred_loglikelihood <- (-2) * (1/n) * (logpi_term + logdensity_term)
  
  return(pred_loglikelihood)
}

load("./Simulations/s1_table.rda")
load("./Simulations/s2_table.rda")
load("./Simulations/s3_table.rda")
