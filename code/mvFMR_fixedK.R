##############################################################################
####              Penalized estimation for a finite mixture of            ####
####                   multivariate regression models                     ####
##############################################################################
####                        Author: Heeyeon Kang                          ####
####                      Supervisor: Sunyoung Shin                       ####
##############################################################################
####         R-code of fitting the data using mvFMR with fixed K          ####
####                        presented in Section 5.                       ####
##############################################################################

# The following R-code fits the simulation data using mvFMR with fixed K
# presented in Section 5.1, 5.2, and 5.3.

### Example ###
# source("./data/simulation_data.R)
# X <- data_generate_5.1.1(12345)$X
# Y <- data_generate_5.1.1(12345)$Y
# mvFMR_fixedK_fit(X, Y, 5)

source("./code/functions.R")

## Assume that K=2 is known ##
mvFMR_fixedK_fit <- function(X, Y, decimal=c("5","6")){
  
  n <- nrow(X)
  P <- ncol(X) - 1
  m <- ncol(Y)
  
  ## Fix K ##
  K <- 2
  
  try({  
    theta_diff <- vector()
    density <- list(list())
    w <- list(list())
    Bk <- list(list())
    sigma <- list(list())
    pi_ <- list(vector())
    
    Bk_array <- array(rep(0, (P+1)*m*K), c(P+1, m, K))
    Bk[[1]] <- lapply(seq(dim(Bk_array)[3]), function(x) Bk_array[ , , x])
    
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
    
    sigma_array <- array(rep(0, m*m*K), c(m, m, K))
    sigma[[1]] <- lapply(lapply(seq(dim(sigma_array)[3]), function(x) sigma_array[ , , x]), function(x) x+diag(m))
    
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
    Bk[[2]] <- m_step_Bk_mvFMR(1)
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
    # Run at least 100 iterations #
    t <- 2
    for(t in 2:100){
      w[[t+1]] <- e_step_w(t)
      
      pi_[[t+1]] <- m_step_pi(t)
      Bk[[t+1]] <- m_step_Bk_mvFMR(t)
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
      
      output <- list(pi=pi_[[t+1]], Bk=Bk[[t+1]], sigma=sigma[[t+1]])
      print(output)
      cat("Iteration :", t+1, "\n")
      
      t <- t+1
    }
    
    while(theta_diff[t] >= 1e-6){
      w[[t+1]] <- e_step_w(t)
      
      pi_[[t+1]] <- m_step_pi(t)
      Bk[[t+1]] <- m_step_Bk_mvFMR(t)
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
      
      output <- list(pi=pi_[[t+1]], Bk=Bk[[t+1]], sigma=sigma[[t+1]])
      print(output)
      cat("Iteration :", t+1, "\n")
      
      if(t == 200){
        break
      }else{
        t <- t+1
      }
    }
  }, silent = TRUE)
  
  optimal <- output
  optimal_w <- w[[t]]
  optimal_density <- density[[t]]
  
  return(list(optimal=optimal, optimal_w=optimal_w, optimal_density=optimal_density))
}
