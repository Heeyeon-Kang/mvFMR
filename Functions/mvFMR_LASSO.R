##############################################################################
####              Penalized estimation for a finite mixture of            ####
####                   multivariate regression models                     ####
##############################################################################
####                        Author: Heeyeon Kang                          ####
####                      Supervisor: Sunyoung Shin                       ####
##############################################################################
####     R-code of fitting the data using mvFMR-LASSO with unknown K      ####
####                      presented in Section 5.3.                       ####
##############################################################################

# The following R-code fits the simulation data using mvFMR-LASSO with unknown K
# presented in Section 5.3.

### Example ###
# source("./Data/simulation_seed_number.R")
# source("./Data/simulation_data.R)
# X <- data_generate_5.3.1(seed_number_5.3.1[1])$X
# Y <- data_generate_5.3.1(seed_number_5.3.1[1])$Y
# mvFMR_LASSO_fit(X, Y, 5)

source("./Functions/functions.R")

## Assume that K is unknown ##
mvFMR_LASSO_fit <- function(X, Y, decimal=c("5","6")){
  
  n <- nrow(X)
  P <- ncol(X) - 1
  m <- ncol(Y)
  
  lambda_s <- c(1e-4, 5*1e-3, 10*1e-3, 15*1e-3, 20*1e-3, 25*1e-3, 
                30*1e-3, 35*1e-3, 40*1e-3, 45*1e-3, 50*1e-3, 55*1e-3, 
                60*1e-3, 65*1e-3, 70*1e-3, 75*1e-3, 80*1e-3, 85*1e-3, 
                90*1e-3, 95*1e-3, 1e-1, 15*1e-2, 2*1e-1, 25*1e-2, 3*1e-1, 
                35*1e-2, 4*1e-1, 45*1e-2, 5*1e-1, 55*1e-2, 6*1e-1, 65*1e-2, 
                7*1e-1, 75*1e-2, 8*1e-1, 85*1e-2, 9*1e-1, 95*1e-2, 1) 
  
  BIC_LASSO_s <- list()
  BIC_LASSO <- vector(length=length(lambda_s))
  total_LASSO <- list(list())
  total_OUTPUTS_LASSO <- list(list())
  optimal_OUTPUTS_LASSO <- list()
  w_s <- list(list())
  optimal_w <- list()
  density_s <- list(list())
  optimal_density <- list()
  
  rho <- 0.4
  e_abs <- 1e-6
  e_rel <- 1e-4
  alpha <- 1.5
  a <- 7
  
  for(K in 1:5){
    for(z in 1:length(lambda_s)){
      
      try({
        theta_diff <- vector()
        density <- list(list())
        w <- list(list())
        Bk <- list(list())
        sigma <- list(list())
        pi_ <- list(vector())
        
        Bk_array <- array(rep(0, (P+1)*m*K), c(P+1, m, K))
        Bk[[1]] <- lapply(seq(dim(Bk_array)[3]), function(x) Bk_array[ , , x])
        
        sigma_array <- array(rep(0, m*m*K), c(m, m, K))
        sigma[[1]] <- lapply(lapply(seq(dim(sigma_array)[3]), function(x) sigma_array[ , , x]), function(x) x+diag(m))
        
        if(K == 1){
          pi_[[1]] <- 1
        }else if(K == 2){
          pi_[[1]] <- c(0.45, 0.55)
        }else if(K == 3){
          pi_[[1]] <- c(0.27, 0.33, 0.4)
        }else if(K == 4){
          pi_[[1]] <- c(0.23, 0.24, 0.26, 0.27)
        }else{
          pi_[[1]] <- c(0.16, 0.18, 0.2, 0.22, 0.24)
        }
        
        lambda <- rep(lambda_s[z], K)
        
        theta_diff[1] <- 0
        
        w_array <- array(rep(0, n*1*K), c(n, 1, K))
        w[[1]] <- lapply(seq(dim(w_array)[3]), function(x) w_array[ , , x])
        
        density[[1]] <- list()
        for(k in 1:K){
          density[[1]][[k]] <- density_f(X,Y,Bk[[1]][[k]],sigma[[1]][[k]])
          for(i in 1:n){
            if(density[[1]][[k]][i] == 0){
              density[[1]][[k]][i] <- 1e-200
            }
          }
        }
        
        ## The first E and M steps ##
        w[[2]] <- e_step_w(1)
        
        pi_[[2]] <- m_step_pi(1)
        Bk[[2]] <- m_step_Bk_mvFMR_L(1)
        sigma[[2]] <- m_step_sigma(1)
        
        density[[2]] <- list()
        for(k in 1:K){
          density[[2]][[k]] <- density_f(X,Y,Bk[[2]][[k]],sigma[[2]][[k]])
          for(i in 1:n){
            if(density[[2]][[k]][i] == 0){
              density[[2]][[k]][i] <- 1e-200
            }
          }
        }
        
        theta_diff[2] <- total_diff_norm(1)
        
        ## The iteration of E and M steps ##
        t <- 2
        while(theta_diff[t] >= 1e-6){
          w[[t+1]] <- e_step_w(t)
          
          pi_[[t+1]] <- m_step_pi(t)
          Bk[[t+1]] <- m_step_Bk_mvFMR_L(t)
          sigma[[t+1]] <- m_step_sigma(t)
          
          density[[t+1]] <- list()
          for(k in 1:K){
            density[[t+1]][[k]] <- density_f(X,Y,Bk[[t+1]][[k]],sigma[[t+1]][[k]])
            for(i in 1:n){
              if(density[[t+1]][[k]][i] == 0){
                density[[t+1]][[k]][i] <- 1e-200
              }
            }
          }
          theta_diff[t+1] <- total_diff_norm(t)
          
          output <- list(lambda=lambda_s[z], pi=pi_[[t+1]], Bk=Bk[[t+1]], sigma=sigma[[t+1]])
          print(output)
          cat("Iteration :", t+1, "\n")
          cat("Lambda :", lambda_s[z], "\n")
          cat("K :", K, "\n")
          
          if(t == 200){
            break
          }else{
            t <- t+1
          }
        }
      }, silent = TRUE)
      
      if(decimal == "5"){
        BIC_LASSO[z] <- BIC_5(t)
      }else if(decimal == "6"){
        BIC_LASSO[z] <- BIC_6(t)
      }else{
        message("Please select 5 or 6")
      }
      
      w_s[[z]] <- w[[t]]
      density_s[[z]] <- density[[t]]
      total_OUTPUTS_LASSO[[z]] <- output
    }
    
    BIC_LASSO_s[[K]] <- replace(BIC_LASSO, which(is.infinite(BIC_LASSO)==T), NA)
    total_LASSO[[K]] <- total_OUTPUTS_LASSO
    optimal_OUTPUTS_LASSO[[K]] <- total_OUTPUTS_LASSO[[which.min(BIC_LASSO_s[[K]])]]
    optimal_w[[K]] <- w_s[[which.min(BIC_LASSO_s[[K]])]]
    optimal_density[[K]] <- density_s[[which.min(BIC_LASSO_s[[K]])]]
  }
  
  optimal_BIC_LASSO_s <- sapply(BIC_LASSO_s, min)
  optimal_K <- which.min(sapply(BIC_LASSO_s, min))
  optimal_OUTPUTS_LASSO_K <- optimal_OUTPUTS_LASSO[[optimal_K]]
  optimal_w_K <- optimal_w[[optimal_K]]
  optimal_density_K <- optimal_density[[optimal_K]]
  
  return(list(total=total_LASSO, optimal=optimal_OUTPUTS_LASSO_K, total_w=optimal_w, 
              total_density=optimal_density, K=optimal_K, optimal_w=optimal_w_K, 
              optimal_density=optimal_density_K, optimal_BIC_LASSO=optimal_BIC_LASSO_s))
}
