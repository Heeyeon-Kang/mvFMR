##############################################################################
####              Penalized estimation for a finite mixture of            ####
####                   multivariate regression models                     ####
##############################################################################
####                        Author: Heeyeon Kang                          ####
####                      Supervisor: Sunyoung Shin                       ####
##############################################################################
####       R-code of fitting the data using the R-package "Flexmix"       ####
####                   presented in Section 5.1 and 5.2.                  ####
##############################################################################

# The following R-code fits the simulation data using the R-package "Flexmix" 
# presented in Section 5.1 and 5.2.

# When fitting using flexmix, there is a problem with cluster shifting 
# compared to the true clusters, making it difficult to calculate 
# the true parameters and MSE.
# Thus, we apply flexmix assuming we know the true parameters, 
# especially, true coefficient matrices, and clusters.

### Example ###
# source("./Data/simulation_seed_number.R")
# source("./Data/simulation_data.R)
# X <- data_generate_5.1.1(seed_number_5.1.1[1])$X
# Y <- data_generate_5.1.1(seed_number_5.1.1[1])$Y
# cluster <- data_generate_5.1.1(seed_number_5.1.1[1])$cluster
# true <- data_generate_5.1.1(seed_number_5.1.1[1])$true
# flexmix_fit(true, cluster, X, Y)

source("./Functions/functions.R")

## Assume that K=2 is known ##
flexmix_fit <- function(true, cluster, X, Y){
  
  n <- nrow(X)
  P <- ncol(X) - 1
  m <- ncol(Y)
  
  ## Fix K ##
  K <- 2
  
  true_Bk <- true$Bk
  total_true_Bk <- matrix(as.vector(sapply(true$Bk, as.vector)), nrow=ncol(X))
  
  apply_flexmix <- list()
  for(i in 1:m){
    apply_flexmix[[i]] <- flexmix(Y[,i] ~ X[,-1], k=2)
  }
  
  clusters_flexmix <- list()
  for(i in 1:m){
    clusters_flexmix[[i]] <- apply_flexmix[[i]]@cluster
  }
  
  rand_index <- vector()
  for(i in 1:m){
    rand_index[i] <- rand.index(cluster, clusters_flexmix[[i]])
  }
  
  cluster_flexmix <- clusters_flexmix[[which.max(rand_index)]]
  
  coef_flexmix <- list()
  for(i in 1:m){
    coef_flexmix[[2*i-1]] <- parameters(apply_flexmix[[i]], component=1, model=1)[1:(P+1)]
    coef_flexmix[[2*i]] <- parameters(apply_flexmix[[i]], component=2, model=1)[1:(P+1)]
  }
  
  total_coef <- sapply(coef_flexmix, cbind)
  
  var_flexmix <- vector()
  for(i in 1:m){
    var_flexmix[2*i-1] <- parameters(apply_flexmix[[i]], component=1, model=1)[P+2]
    var_flexmix[2*i] <- parameters(apply_flexmix[[i]], component=2, model=1)[P+2]
  }
  total_var <- var_flexmix
  
  norm_of_true_and_pred_Bk <- vector()
  
  list_index <- as.list(as.data.frame(matrix(1:(m*K), nrow=K)))
  index <- as.matrix(expand.grid(list_index))
  
  total_coef <- list()
  for(j in 1:nrow(index)){
    total_coef[[j]] <- list(matrix(nrow=(P+1), ncol=m), matrix(nrow=(P+1), ncol=m))
    total <- 1:(m*K)
    id1 <- as.vector(index[j,])
    id2 <- total[-id1]
    for(i in 1:m){
      total_coef[[j]][[1]][,i] <- coef_flexmix[[id1[i]]]
      total_coef[[j]][[2]][,i] <- coef_flexmix[[id2[i]]]
    }
  }
  
  total_big_coef <- list()
  total_big_coef <- lapply(total_coef, function(x) do.call(cbind, x))
  
  for(j in 1:nrow(index)){
    norm_of_true_and_pred_Bk[j] <- norm(total_true_Bk - total_big_coef[[j]], "F")
  }
  
  switch_of_pred <- as.vector(index[which.min(norm_of_true_and_pred_Bk),])
  other_switch_of_pred <- total[-switch_of_pred]
  
  switched_total_coef <- total_coef[[which.min(norm_of_true_and_pred_Bk)]]
  switched_total_var <- c(total_var[switch_of_pred], total_var[other_switch_of_pred])
  
  flexmix_pi <- as.vector(table(cluster_flexmix)) / n
  flexmix_Bk <- switched_total_coef
  flexmix_sigma <- list(diag(total_var[1:m]), diag(total_var[-(1:m)]))
  
  pi_diff <- vector()
  for(i in 1:nrow(permutations(K,K))){
    pi_diff[i] <- norm(true$pi - flexmix_pi[permutations(K,K)[i,]], "2")
  }
  MSE_pi <- 1/K * (min(pi_diff)^2)
  
  MSE_Bk <- 1/K * (norm(total_true_Bk - total_big_coef[[which.min(norm_of_true_and_pred_Bk)]], "F")^2)
  
  total_pred_sigma <- sapply(flexmix_sigma, as.vector)
  total_true_sigma <- sapply(true$sigma, as.vector)
  MSE_sigma <- 1/K * (norm(total_pred_sigma - total_true_sigma, "F")^2)
  output <- list(pi=flexmix_pi, Bk=flexmix_Bk, sigma=flexmix_sigma)
  
  return(list(optimal=output, cluster=cluster_flexmix, MSE=c(MSE_pi, MSE_Bk, MSE_sigma)))
}
