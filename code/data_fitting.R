##############################################################################
####              Penalized estimation for a finite mixture of            ####
####                   multivariate regression models                     ####
##############################################################################
####                        Author: Heeyeon Kang                          ####
####                      Supervisor: Sunyoung Shin                       ####
##############################################################################
####           R-code of fitting the data presented in Section 5.         ####
##############################################################################

# The following R-code contains all functions of fitting the simulation data
# presented in Section 5.1, 5.2, and 5.3.

# The code consists of functions for fitting the data using our method 
# when K is known and not, fitting with the R-package flexmix, 
# and obtaining the oracle estimator.

# When fitting using flexmix, there is a problem with cluster shifting 
# compared to the true clusters, making it difficult to calculate 
# the true parameters and MSE.
# Thus, we apply flexmix assuming we know the true parameters, 
# especially, true coefficient matrices, and clusters.

source("./code/functions.R")

### Fit the data using our method, mvFMR ###
## Assume that K=2 is known with penalty function ##
# mvFMR with no penalty function #
mvFMR_givenK_fit <- function(X, Y, decimal=c("5","6")){
  
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

# mvFMR with LASSO penalty function #
mvFMR_LASSO_givenK_fit <- function(X, Y, decimal=c("5","6")){

  n <- nrow(X)
  P <- ncol(X) - 1
  m <- ncol(Y)
  
  ## Fix K ##
  K <- 2
  
  lambda_s <- c(1e-4, 5*1e-3, 10*1e-3, 15*1e-3, 20*1e-3, 25*1e-3, 
                30*1e-3, 35*1e-3, 40*1e-3, 45*1e-3, 50*1e-3, 55*1e-3, 
                60*1e-3, 65*1e-3, 70*1e-3, 75*1e-3, 80*1e-3, 85*1e-3, 
                90*1e-3, 95*1e-3, 1e-1, 15*1e-2, 2*1e-1, 25*1e-2, 3*1e-1, 
                35*1e-2, 4*1e-1, 45*1e-2, 5*1e-1, 55*1e-2, 6*1e-1, 65*1e-2, 
                7*1e-1, 75*1e-2, 8*1e-1, 85*1e-2, 9*1e-1, 95*1e-2, 1) 
  
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
      message("Please select 5 or 6.")
    }
    
    w_s[[z]] <- w[[t]]
    density_s[[z]] <- density[[t]]
    total_OUTPUTS_LASSO[[z]] <- output
  }
  
  optimal_OUTPUTS_LASSO <- total_OUTPUTS_LASSO[[which.min(BIC_LASSO)]]
  optimal_w <- w_s[[which.min(BIC_LASSO)]]
  optimal_density <- density_s[[which.min(BIC_LASSO)]]
  
  return(list(total=total_OUTPUTS_LASSO, optimal=optimal_OUTPUTS_LASSO, total_w=w_s, 
                 total_density=density_s, optimal_w=optimal_w, optimal_density=optimal_density))
}

# mvFMR with SCAD penalty function #
mvFMR_SCAD_givenK_fit <- function(X, Y, decimal=c("5","6")){
  
  n <- nrow(X)
  P <- ncol(X) - 1
  m <- ncol(Y)
  K <- 2 # Fix K
  
  lambda_s <- c(1e-4, 5*1e-3, 10*1e-3, 15*1e-3, 20*1e-3, 25*1e-3, 
                30*1e-3, 35*1e-3, 40*1e-3, 45*1e-3, 50*1e-3, 55*1e-3, 
                60*1e-3, 65*1e-3, 70*1e-3, 75*1e-3, 80*1e-3, 85*1e-3, 
                90*1e-3, 95*1e-3, 1e-1, 15*1e-2, 2*1e-1, 25*1e-2, 3*1e-1, 
                35*1e-2, 4*1e-1, 45*1e-2, 5*1e-1, 55*1e-2, 6*1e-1, 65*1e-2, 
                7*1e-1, 75*1e-2, 8*1e-1, 85*1e-2, 9*1e-1, 95*1e-2, 1) 
  
  BIC_SCAD <- vector(length=length(lambda_s))
  total_OUTPUTS_SCAD <- list(list())
  optimal_OUTPUTS_SCAD <- list()
  w_s <- list(list())
  optimal_w <- list()
  density_s <- list(list())
  optimal_density <- list()
  
  rho <- 0.4
  e_abs <- 1e-6
  e_rel <- 1e-4
  alpha <- 1.5
  a <- 7
  
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
      Bk[[2]] <- m_step_Bk_mvFMR_S(1)
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
        Bk[[t+1]] <- m_step_Bk_mvFMR_S(t)
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
      BIC_SCAD[z] <- BIC_5(t)
    }else if(decimal == "6"){
      BIC_SCAD[z] <- BIC_6(t)
    }else{
      message("Please select 5 or 6")
    }
    
    total_OUTPUTS_SCAD[[z]] <- output
    w_s[[z]] <- w[[t]]
    density_s[[z]] <- density[[t]]
  }
  
  optimal_OUTPUTS_SCAD <- total_OUTPUTS_SCAD[[which.min(BIC_SCAD)]]
  optimal_w <- w_s[[which.min(BIC_SCAD)]]
  optimal_density <- density_s[[which.min(BIC_SCAD)]]
  
  return(list(total=total_OUTPUTS_SCAD, optimal=optimal_OUTPUTS_SCAD, total_w=w_s, 
              total_density=density_s, optimal_w=optimal_w, 
              optimal_density=optimal_density))
}

# mvFMR with MCP penalty function #
mvFMR_MCP_givenK_fit <- function(X, Y, decimal=c("5","6")){
  
  n <- nrow(X)
  P <- ncol(X) - 1
  m <- ncol(Y)
  K <- 2 # Fix K
  
  lambda_s <- c(1e-4, 5*1e-3, 10*1e-3, 15*1e-3, 20*1e-3, 25*1e-3, 
                30*1e-3, 35*1e-3, 40*1e-3, 45*1e-3, 50*1e-3, 55*1e-3, 
                60*1e-3, 65*1e-3, 70*1e-3, 75*1e-3, 80*1e-3, 85*1e-3, 
                90*1e-3, 95*1e-3, 1e-1, 15*1e-2, 2*1e-1, 25*1e-2, 3*1e-1, 
                35*1e-2, 4*1e-1, 45*1e-2, 5*1e-1, 55*1e-2, 6*1e-1, 65*1e-2, 
                7*1e-1, 75*1e-2, 8*1e-1, 85*1e-2, 9*1e-1, 95*1e-2, 1) 
  
  BIC_MCP <- vector(length=length(lambda_s))
  total_OUTPUTS_MCP <- list(list())
  optimal_OUTPUTS_MCP <- list()
  w_s <- list(list())
  optimal_w <- list()
  density_s <- list(list())
  optimal_density <- list()
  
  rho <- 0.4
  e_abs <- 1e-6
  e_rel <- 1e-4
  alpha <- 1.5
  a <- 7
  
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
      Bk[[2]] <- m_step_Bk_mvFMR_M(1)
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
        Bk[[t+1]] <- m_step_Bk_mvFMR_M(t)
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
      BIC_MCP[z] <- BIC_5(t)
    }else if(decimal == "6"){
      BIC_MCP[z] <- BIC_6(t)
    }else{
      message("Please select 5 or 6")
    }
    
    w_s[[z]] <- w[[t]]
    density_s[[z]] <- density[[t]]
    total_OUTPUTS_MCP[[z]] <- output
  }
  
  optimal_OUTPUTS_MCP <- total_OUTPUTS_MCP[[which.min(BIC_MCP)]]
  optimal_w <- w_s[[which.min(BIC_MCP)]]
  optimal_density <- density_s[[which.min(BIC_MCP)]]
  
  return(list(total=total_OUTPUTS_MCP, optimal=optimal_OUTPUTS_MCP, total_w=w_s,
              total_density=density_s, optimal_w=optimal_w, 
              optimal_density=optimal_density))
}


## Assume that K is unknown with penalty function ##
# mvFMR with no penalty function #
mvFMR_fit <- function(X, Y, decimal=c("5","6")){
  
  n <- nrow(X)
  P <- ncol(X) - 1
  m <- ncol(Y)
  
  w_s <- list()
  density_s <- list()
  total_OUTPUTS_UNPEN <- list()
  BIC_UNPEN <- vector()
  
  for(K in 1:5){
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
        cat("K :", K, "\n")
        
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
        cat("K :", K, "\n")
        
        if(t == 200){
          break
        }else{
          t <- t+1
        }
      }
    }, silent = TRUE)
    
    total_OUTPUTS_UNPEN[[K]] <- output
    BIC_UNPEN[K] <- BIC_5(t)
    w_s[[K]] <- w[[t]]
    density_s[[K]] <- density[[t]]
  }
  
  BIC_UNPEN <- replace(BIC_UNPEN, which(is.infinite(BIC_UNPEN)==T), NA)
  optimal_K <- which.min(BIC_UNPEN)
  optimal_OUTPUTS_UNPEN <- total_OUTPUTS_UNPEN[[optimal_K]]
  optimal_w <- w_s[[optimal_K]]
  optimal_denstiy <- density_s[[optimal_K]]
  
  return(list(optimal=optimal_OUTPUTS_UNPEN, K=optimal_K, total_w=w_s, total_density=density_s,
              optimal_w=optimal_w, optimal_density=optimal_density))
}

# mvFMR with LASSO penalty function #
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

# mvFMR with SCAD penalty function #
mvFMR_SCAD_fit <- function(X, Y, decimal=c("5","6")){
  
  n <- nrow(X)
  P <- ncol(X) - 1
  m <- ncol(Y)
  
  lambda_s <- c(1e-4, 5*1e-3, 10*1e-3, 15*1e-3, 20*1e-3, 25*1e-3, 
                30*1e-3, 35*1e-3, 40*1e-3, 45*1e-3, 50*1e-3, 55*1e-3, 
                60*1e-3, 65*1e-3, 70*1e-3, 75*1e-3, 80*1e-3, 85*1e-3, 
                90*1e-3, 95*1e-3, 1e-1, 15*1e-2, 2*1e-1, 25*1e-2, 3*1e-1, 
                35*1e-2, 4*1e-1, 45*1e-2, 5*1e-1, 55*1e-2, 6*1e-1, 65*1e-2, 
                7*1e-1, 75*1e-2, 8*1e-1, 85*1e-2, 9*1e-1, 95*1e-2, 1) 
  
  BIC_SCAD_s <- list()
  BIC_SCAD <- vector(length=length(lambda_s))
  total_SCAD <- list(list())
  total_OUTPUTS_SCAD <- list(list())
  optimal_OUTPUTS_SCAD <- list()
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
        Bk[[2]] <- m_step_Bk_mvFMR_S(1)
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
          Bk[[t+1]] <- m_step_Bk_mvFMR_S(t)
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
        BIC_SCAD[z] <- BIC_5(t)
      }else if(decimal == "6"){
        BIC_SCAD[z] <- BIC_6(t)
      }else{
        message("Please select 5 or 6")
      }
      
      total_OUTPUTS_SCAD[[z]] <- output
      w_s[[z]] <- w[[t]]
      density_s[[z]] <- density[[t]]
    }
    
    BIC_SCAD_s[[K]] <- replace(BIC_SCAD, which(is.infinite(BIC_SCAD)==T), NA)
    total_SCAD[[K]] <- total_OUTPUTS_SCAD
    optimal_OUTPUTS_SCAD[[K]] <- total_OUTPUTS_SCAD[[which.min(BIC_SCAD_s[[K]])]]
    optimal_w[[K]] <- w_s[[which.min(BIC_SCAD_s[[K]])]]
    optimal_density[[K]] <- density_s[[which.min(BIC_SCAD_s[[K]])]]
  }
  
  optimal_BIC_SCAD_s <- sapply(BIC_SCAD_s, min)
  optimal_K <- which.min(sapply(BIC_SCAD_s, min))
  optimal_OUTPUTS_SCAD_K <- optimal_OUTPUTS_SCAD[[optimal_K]]
  optimal_w_K <- optimal_w[[optimal_K]]
  optimal_density_K <- optimal_density[[optimal_K]]
  
  return(list(total=total_SCAD, optimal=optimal_OUTPUTS_SCAD_K, total_w=optimal_w, 
              total_density=optimal_density, K=optimal_K, optimal_w=optimal_w_K, 
              optimal_density=optimal_density_K, optimal_BIC_SCAD=optimal_BIC_SCAD_s))
}

# mvFMR with MCP penalty function #
mvFMR_MCP_fit <- function(X, Y, decimal=c("5","6")){
  
  n <- nrow(X)
  P <- ncol(X) - 1
  m <- ncol(Y)
  
  lambda_s <- c(1e-4, 5*1e-3, 10*1e-3, 15*1e-3, 20*1e-3, 25*1e-3, 
                30*1e-3, 35*1e-3, 40*1e-3, 45*1e-3, 50*1e-3, 55*1e-3, 
                60*1e-3, 65*1e-3, 70*1e-3, 75*1e-3, 80*1e-3, 85*1e-3, 
                90*1e-3, 95*1e-3, 1e-1, 15*1e-2, 2*1e-1, 25*1e-2, 3*1e-1, 
                35*1e-2, 4*1e-1, 45*1e-2, 5*1e-1, 55*1e-2, 6*1e-1, 65*1e-2, 
                7*1e-1, 75*1e-2, 8*1e-1, 85*1e-2, 9*1e-1, 95*1e-2, 1) 
  
  BIC_MCP_s <- list()
  BIC_MCP <- vector(length=length(lambda_s))
  total_MCP <- list(list())
  total_OUTPUTS_MCP <- list(list())
  optimal_OUTPUTS_MCP <- list()
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
        Bk[[2]] <- m_step_Bk_mvFMR_M(1)
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
          Bk[[t+1]] <- m_step_Bk_mvFMR_M(t)
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
        BIC_MCP[z] <- BIC_5(t)
      }else if(decimal == "6"){
        BIC_MCP[z] <- BIC_6(t)
      }else{
        message("Please select 5 or 6")
      }
      
      w_s[[z]] <- w[[t]]
      density_s[[z]] <- density[[t]]
      total_OUTPUTS_MCP[[z]] <- output
    }
    
    BIC_MCP_s[[K]] <- replace(BIC_MCP, which(is.infinite(BIC_MCP)==T), NA)
    total_MCP[[K]] <- total_OUTPUTS_MCP
    optimal_OUTPUTS_MCP[[K]] <- total_OUTPUTS_MCP[[which.min(BIC_MCP_s[[K]])]]
    optimal_w[[K]] <- w_s[[which.min(BIC_MCP_s[[K]])]]
    optimal_density[[K]] <- density_s[[which.min(BIC_MCP_s[[K]])]]
  }
  
  optimal_BIC_MCP_s <- sapply(BIC_MCP_s, min)
  optimal_K <- which.min(sapply(BIC_MCP_s, min))
  optimal_OUTPUTS_MCP_K <- optimal_OUTPUTS_MCP[[optimal_K]]
  optimal_w_K <- optimal_w[[optimal_K]]
  optimal_density_K <- optimal_density[[optimal_K]]
  
  return(list(total=total_MCP, optimal=optimal_OUTPUTS_MCP_K, total_w=optimal_w, 
              total_density=optimal_density, K=optimal_K, optimal_w=optimal_w_K, 
              optimal_density=optimal_density_K, optimal_BIC_MCP=optimal_BIC_MCP_s))
}


### Fit the data using the R-package, flexmix ###
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
    pi_diff[i] <- norm(true_pi_ - flexmix_pi[permutations(K,K)[i,]], "2")
  }
  MSE_pi <- 1/K * (min(pi_diff)^2)
  
  MSE_Bk <- 1/K * (norm(total_true_Bk - total_big_coef[[which.min(norm_of_true_and_pred_Bk)]], "F")^2)
  
  total_pred_sigma <- sapply(flexmix_sigma, as.vector)
  total_true_sigma <- sapply(true$sigma, as.vector)
  MSE_sigma <- 1/K * (norm(total_pred_sigma - total_true_sigma, "F")^2)
  output <- list(pi=flexmix_pi, Bk=flexmix_Bk, sigma=flexmix_sigma)
  
  return(list(optimal=output, cluster=cluster_flexmix, MSE=c(MSE_pi, MSE_Bk, MSE_sigma)))
}

### Obtain the oracle estimator ###
oracle_fit <- function(X, Y){
  
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
          density[[1]][[k]][i] <- 1e-323
        }
      }
    }
    
    ## The first E and M steps ##
    w[[2]] <- e_step_w(1)
    
    pi_[[2]] <- m_step_pi(1)
    Bk[[2]] <- m_step_Bk_oracle(1)
    sigma[[2]] <- m_step_sigma(1)
    
    density[[2]] <- list()
    for(k in 1:K){
      density[[2]][[k]] <- density_f(X,Y,Bk[[2]][[k]],sigma[[2]][[k]])
      for(i in 1:n){
        if(density[[2]][[k]][i] == 0){
          density[[2]][[k]][i] <- 1e-323
        }
      }
    }
    
    theta_diff[2] <- total_diff_norm(1)
    
    ## The iteration of E and M steps ##
    t <- 2
    for(t in 2:50){
      w[[t+1]] <- e_step_w(t)
      
      pi_[[t+1]] <- m_step_pi(t)
      Bk[[t+1]] <- m_step_Bk_oracle(t)
      sigma[[t+1]] <- m_step_sigma(t)
      
      density[[t+1]] <- list()
      for(k in 1:K){
        density[[t+1]][[k]] <- density_f(X,Y,Bk[[t+1]][[k]],sigma[[t+1]][[k]])
        for(i in 1:n){
          if(density[[t+1]][[k]][i] == 0){
            density[[t+1]][[k]][i] <- 1e-323
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
      Bk[[t+1]] <- m_step_Bk_oracle(t)
      sigma[[t+1]] <- m_step_sigma(t)
      
      density[[t+1]] <- list()
      for(k in 1:K){
        density[[t+1]][[k]] <- density_f(X,Y,Bk[[t+1]][[k]],sigma[[t+1]][[k]])
        for(i in 1:n){
          if(density[[t+1]][[k]][i] == 0){
            density[[t+1]][[k]][i] <- 1e-323
          }
        }
      }
      theta_diff[t+1] <- total_diff_norm(t)
      
      output <- list(pi=pi_[[t+1]], Bk=Bk[[t+1]], sigma=sigma[[t+1]])
      print(output)
      cat("Iteration :", t+1, "\n")
      
      if(t == 100){
        break
      }else{
        t <- t+1
      }
    }
  }, silent = TRUE)
  
  return(list(optimal=output, w=w[[t]], density=density[[t]]))
}

