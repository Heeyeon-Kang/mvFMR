##############################################################################
####              Penalized estimation for a finite mixture of            ####
####                   multivariate regression models                     ####
##############################################################################
####                        Author: Heeyeon Kang                          ####
####                      Supervisor: Sunyoung Shin                       ####
##############################################################################
####           R-code of functions used for EM-ADMM algorithm.            ####
##############################################################################

# The following R-code contains all functions for running EM-ADMM algorithm.


### Packages install ###
requiredPackages <- c("MASS", "pracma", "expm", "corrcoverage", 
                      "causact", "fossil", "flexmix", "gtools")
for(p in requiredPackages){
  if(!require(p,character.only = TRUE)) install.packages(p)
  library(p,character.only = TRUE)
}

### The density of multivariate Gaussian regression model ###
density_f <- function(X,Y,B,sigma){
  
  n <- nrow(X)
  m <- ncol(Y)
  
  a <- (2*pi)^(-m/2)
  b <- determinant(sigma,logarithm=FALSE)$sign
  c <- determinant(sigma,logarithm=FALSE)$modulus[1]
  d <- (b*c)^(-1/2)
  e <- vector(length=n)
  f <- vector(length=n)
  for(i in 1:n){
    e[i] <- (-1/2)*(t(Y[i,]-t(B)%*%X[i,])%*%ginv(sigma)%*%(Y[i,]-t(B)%*%X[i,]))
    f[i] <- a*d*exp(e[i])
  }
  return(f)
}

### ADMM functions by penalty type ###
admm_Bk_LASSO <- function(l){
  
  Bk_sol <- list(list()) # Solution of Bk after applying ADMM.
  Ck <- list(list())
  Hk <- list(list())
  U1 <- list(list())
  epsilon <- list(list())
  R <- list(list()) # Primal residual
  S <- list(list()) # Dual residual
  
  ## Initialization ##
  Bk_sol[[1]] <- Bk[[l]]
  Ck[[1]] <- list()
  U1[[1]] <- list()
  Hk[[1]] <- list()
  for(k in 1:K){
    Ck[[1]][[k]] <- Bk_sol[[1]][[k]]
    U1[[1]][[k]] <- matrix(0, nrow=P+1, ncol=m)
  }
  
  ## Sylvester's equation ##
  A_array <- array(rep(0, (P+1)*(P+1)*K), c(P+1, P+1, K))
  A <- lapply(seq(dim(A_array)[3]), function(x) A_array[ , , x])
  for(k in 1:K){
    for(i in 1:n){
      A[[k]] <- A[[k]] + (1/n)*w[[l+1]][[k]][i]*(X[i,]%*%t(X[i,]))
    }
  }
  
  E <- list()
  for(k in 1:K){
    E[[k]] <- -rho * sigma[[l]][[k]]
  }
  
  W_array <- array(rep(0, (P+1)*m*K), c(P+1, m, K))
  W <- lapply(seq(dim(W_array)[3]), function(x) W_array[ , , x]) 
  for(k in 1:K){
    for(i in 1:n){
      W[[k]] <- W[[k]] + (1/n)*w[[l+1]][[k]][i]*(X[i,]%*%t(Y[i,]))
    }
  }
  
  G <- list()
  for(k in 1:K){
    G[[k]] <- W[[k]] + rho*(Ck[[1]][[k]]-U1[[1]][[k]])%*%sigma[[l]][[k]]
  }
  
  delta <- list()
  for(k in 1:K){
    delta[[k]] <- charpoly(A[[k]])
  }
  
  fAE_array <- array(rep(0, m*m*K), c(m, m, K))
  fAE <- lapply(seq(dim(fAE_array)[3]), function(x) fAE_array[ , , x])
  for(k in 1:K){
    for(j in 1:(P+2)){
      fAE[[k]] <- fAE[[k]] + delta[[k]][j]*(E[[k]]%^%(P+2-j))
    }
  }
  
  sumAGE <- as.list(array(rep(0, K), c(1, 1, K)))
  for(k in 1:K){
    for(p in 1:(P+1)){
      for(q in 0:(p-1)){
        if(p-q-1 >= 0){
          sumAGE[[k]] <- sumAGE[[k]] + delta[[k]][P+2-p]*((A[[k]]%^%(q))%*%G[[k]]%*%(E[[k]]%^%(p-q-1)))
        }
      }
    }
  }
  
  ## Update Bk ##
  Bk_sol[[2]] <- list()
  for(k in 1:K){
    Bk_sol[[2]][[k]] <- -sumAGE[[k]]%*%ginv(fAE[[k]])
  }
  
  ## Update Ck ##
  Ck[[2]] <- list()
  Hk[[2]] <- list()
  for(k in 1:K){
    Hk[[2]][[k]] <- alpha*Bk_sol[[2]][[k]] + (1-alpha)*Ck[[1]][[k]]
  }
  for(k in 1:K){
    Ck[[2]][[k]] <- LASSO_thresholding(Hk[[2]][[k]]+U1[[1]][[k]], pi_[[l+1]][k], lambda[k], rho)
  }
  
  ## Update U1 ##
  U1[[2]] <- list()
  for(k in 1:K){
    U1[[2]][[k]] <- U1[[1]][[k]] + Hk[[2]][[k]] - Ck[[2]][[k]]
  } 
  
  ## Residuals ##
  R[[2]] <- list()
  S[[2]] <- list()
  for(k in 1:K){
    R[[2]][[k]] <- primal_residual(Bk_sol[[2]][[k]], Ck[[2]][[k]])
    S[[2]][[k]] <- dual_residual(U1[[2]][[k]], U1[[1]][[k]], rho)
  }
  
  epsilon[[2]] <- matrix(nrow=K, ncol=2)
  for(k in 1:K){
    e_pri <- epsilon_primal(Bk_sol[[2]][[k]], Ck[[2]][[k]], e_abs,e_rel)
    e_dual <- epsilon_dual(U1[[2]][[k]], e_abs, e_rel)
    epsilon[[2]][k,] <- c(e_pri, e_dual)
  }
  
  ## The iterations of updating Bk, Ck, U1 ##
  t <- 2
  while(stop_admm_lasso(R[[t]], S[[t]], epsilon[[t]])){
    
    for(k in 1:K){
      G[[k]] <- W[[k]] + rho*(Ck[[t]][[k]]-U1[[t]][[k]])%*%sigma[[l]][[k]]
    }
    
    sumAGE <- as.list(array(rep(0, K), c(1, 1, K)))
    for(k in 1:K){
      for(p in 1:(P+1)){
        for(q in 0:(p-1)){
          if(p-q-1 >= 0){
            sumAGE[[k]] <- sumAGE[[k]] + 
              delta[[k]][P+2-p]*((A[[k]]%^%(q))%*%G[[k]]%*%(E[[k]]%^%(p-q-1)))
          }
        }
      }
    }
    
    ## Update Bk ##
    Bk_sol[[t+1]] <- list()
    for(k in 1:K){
      Bk_sol[[t+1]][[k]] <- -sumAGE[[k]]%*%ginv(fAE[[k]])
    }
    
    ## Update Ck ##
    Ck[[t+1]] <- list()
    Hk[[t+1]] <- list()
    for(k in 1:K){
      Hk[[t+1]][[k]] <- alpha*Bk_sol[[t+1]][[k]] + (1-alpha)*Ck[[t]][[k]]
    }
    for(k in 1:K){
      Ck[[t+1]][[k]] <- LASSO_thresholding(Hk[[t+1]][[k]]+U1[[t]][[k]], pi_[[l+1]][k], lambda[k], rho)
    }
    
    ## Update U1 ##
    U1[[t+1]] <- list()
    for(k in 1:K){
      U1[[t+1]][[k]] <- U1[[t]][[k]] + Hk[[t+1]][[k]] - Ck[[t+1]][[k]]
    } 
    
    ## Residuals ##
    R[[t+1]] <- list()
    S[[t+1]] <- list()
    for(k in 1:K){
      R[[t+1]][[k]] <- primal_residual(Bk_sol[[t+1]][[k]], Ck[[t+1]][[k]])
      S[[t+1]][[k]] <- dual_residual(U1[[t+1]][[k]], U1[[t]][[k]], rho)
    }
    
    epsilon[[t+1]] <- matrix(nrow=K, ncol=2)
    for(k in 1:K){
      e_pri <- epsilon_primal(Bk_sol[[t+1]][[k]], Ck[[t+1]][[k]], e_abs, e_rel)
      e_dual <- epsilon_dual(U1[[t+1]][[k]], e_abs, e_rel)
      epsilon[[t+1]][k,] <- c(e_pri, e_dual)
    }
    
    if(t == 10000){
      break
    }else{
      t <- t+1
    }
  }
  return(Bk_sol[[t]])
}

admm_Bk_SCAD <- function(l){
  
  Bk_sol <- list(list()) # Solution of Bk after applying ADMM.
  Ck <- list(list())
  Hk <- list(list())
  U1 <- list(list())
  epsilon <- list(list())
  R <- list(list()) # Primal residual
  S <- list(list()) # Dual residual
  
  ## Initialization ##
  Bk_sol[[1]] <- Bk[[l]]
  Ck[[1]] <- list()
  U1[[1]] <- list()
  Hk[[1]] <- list()
  for(k in 1:K){
    Ck[[1]][[k]] <- Bk_sol[[1]][[k]]
    U1[[1]][[k]] <- matrix(0, nrow=P+1, ncol=m)
  }
  
  ## Sylvester's equation ##
  A_array <- array(rep(0, (P+1)*(P+1)*K), c(P+1, P+1, K))
  A <- lapply(seq(dim(A_array)[3]), function(x) A_array[ , , x])
  for(k in 1:K){
    for(i in 1:n){
      A[[k]] <- A[[k]] + (1/n)*w[[l+1]][[k]][i]*(X[i,]%*%t(X[i,]))
    }
  }
  
  E <- list()
  for(k in 1:K){
    E[[k]] <- -rho * sigma[[l]][[k]]
  }
  
  W_array <- array(rep(0, (P+1)*m*K), c(P+1, m, K))
  W <- lapply(seq(dim(W_array)[3]), function(x) W_array[ , , x]) 
  for(k in 1:K){
    for(i in 1:n){
      W[[k]] <- W[[k]] + (1/n)*w[[l+1]][[k]][i]*(X[i,]%*%t(Y[i,]))
    }
  }
  
  G <- list()
  for(k in 1:K){
    G[[k]] <- W[[k]] + rho*(Ck[[1]][[k]]-U1[[1]][[k]])%*%sigma[[l]][[k]]
  }
  
  delta <- list()
  for(k in 1:K){
    delta[[k]] <- charpoly(A[[k]])
  }
  
  fAE_array <- array(rep(0, m*m*K), c(m, m, K))
  fAE <- lapply(seq(dim(fAE_array)[3]), function(x) fAE_array[ , , x])
  for(k in 1:K){
    for(j in 1:(P+2)){
      fAE[[k]] <- fAE[[k]] + delta[[k]][j]*(E[[k]]%^%(P+2-j))
    }
  }
  
  sumAGE <- as.list(array(rep(0, K), c(1, 1, K)))
  for(k in 1:K){
    for(p in 1:(P+1)){
      for(q in 0:(p-1)){
        if(p-q-1 >= 0){
          sumAGE[[k]] <- sumAGE[[k]] + delta[[k]][P+2-p]*((A[[k]]%^%(q))%*%G[[k]]%*%(E[[k]]%^%(p-q-1)))
        }
      }
    }
  }
  
  ## Update Bk ##
  Bk_sol[[2]] <- list()
  for(k in 1:K){
    Bk_sol[[2]][[k]] <- -sumAGE[[k]]%*%ginv(fAE[[k]])
  }
  
  ## Update Ck ##
  Ck[[2]] <- list()
  Hk[[2]] <- list()
  for(k in 1:K){
    Hk[[2]][[k]] <- alpha*Bk_sol[[2]][[k]] + (1-alpha)*Ck[[1]][[k]]
  }
  for(k in 1:K){
    Ck[[2]][[k]] <- SCAD_thresholding(Hk[[2]][[k]]+U1[[1]][[k]], pi_[[l+1]][k], lambda[k], rho, a)
  }
  
  ## Update U1 ##
  U1[[2]] <- list()
  for(k in 1:K){
    U1[[2]][[k]] <- U1[[1]][[k]] + Hk[[2]][[k]] - Ck[[2]][[k]]
  } 
  
  ## Residuals ##
  R[[2]] <- list()
  S[[2]] <- list()
  for(k in 1:K){
    R[[2]][[k]] <- primal_residual(Bk_sol[[2]][[k]], Ck[[2]][[k]])
    S[[2]][[k]] <- dual_residual(U1[[2]][[k]], U1[[1]][[k]], rho)
  }
  
  epsilon[[2]] <- matrix(nrow=K, ncol=2)
  for(k in 1:K){
    e_pri <- epsilon_primal(Bk_sol[[2]][[k]], Ck[[2]][[k]], e_abs,e_rel)
    e_dual <- epsilon_dual(U1[[2]][[k]], e_abs, e_rel)
    epsilon[[2]][k,] <- c(e_pri, e_dual)
  }
  
  ## The iterations of updating Bk, Ck, U1 ##
  t <- 2
  while(stop_admm_lasso(R[[t]], S[[t]], epsilon[[t]])){
    
    for(k in 1:K){
      G[[k]] <- W[[k]] + rho*(Ck[[t]][[k]]-U1[[t]][[k]])%*%sigma[[l]][[k]]
    }
    
    sumAGE <- as.list(array(rep(0, K), c(1, 1, K)))
    for(k in 1:K){
      for(p in 1:(P+1)){
        for(q in 0:(p-1)){
          if(p-q-1 >= 0){
            sumAGE[[k]] <- sumAGE[[k]] + 
              delta[[k]][P+2-p]*((A[[k]]%^%(q))%*%G[[k]]%*%(E[[k]]%^%(p-q-1)))
          }
        }
      }
    }
    
    ## Update Bk ##
    Bk_sol[[t+1]] <- list()
    for(k in 1:K){
      Bk_sol[[t+1]][[k]] <- -sumAGE[[k]]%*%ginv(fAE[[k]])
    }
    
    ## Update Ck ##
    Ck[[t+1]] <- list()
    Hk[[t+1]] <- list()
    for(k in 1:K){
      Hk[[t+1]][[k]] <- alpha*Bk_sol[[t+1]][[k]] + (1-alpha)*Ck[[t]][[k]]
    }
    for(k in 1:K){
      Ck[[t+1]][[k]] <- SCAD_thresholding(Hk[[t+1]][[k]]+U1[[t]][[k]], pi_[[l+1]][k], lambda[k], rho, a)
    }
    
    ## Update U1 ##
    U1[[t+1]] <- list()
    for(k in 1:K){
      U1[[t+1]][[k]] <- U1[[t]][[k]] + Hk[[t+1]][[k]] - Ck[[t+1]][[k]]
    } 
    
    ## Residuals ##
    R[[t+1]] <- list()
    S[[t+1]] <- list()
    for(k in 1:K){
      R[[t+1]][[k]] <- primal_residual(Bk_sol[[t+1]][[k]], Ck[[t+1]][[k]])
      S[[t+1]][[k]] <- dual_residual(U1[[t+1]][[k]], U1[[t]][[k]], rho)
    }
    
    epsilon[[t+1]] <- matrix(nrow=K, ncol=2)
    for(k in 1:K){
      e_pri <- epsilon_primal(Bk_sol[[t+1]][[k]], Ck[[t+1]][[k]], e_abs, e_rel)
      e_dual <- epsilon_dual(U1[[t+1]][[k]], e_abs, e_rel)
      epsilon[[t+1]][k,] <- c(e_pri, e_dual)
    }
    
    if(t == 10000){
      break
    }else{
      t <- t+1
    }
  }
  return(Bk_sol[[t]])
}

admm_Bk_MCP <- function(l){
  
  Bk_sol <- list(list()) # Solution of Bk after applying ADMM.
  Ck <- list(list())
  Hk <- list(list())
  U1 <- list(list())
  epsilon <- list(list())
  R <- list(list()) # Primal residual
  S <- list(list()) # Dual residual
  
  ## Initialization ##
  Bk_sol[[1]] <- Bk[[l]]
  Ck[[1]] <- list()
  U1[[1]] <- list()
  Hk[[1]] <- list()
  for(k in 1:K){
    Ck[[1]][[k]] <- Bk_sol[[1]][[k]]
    U1[[1]][[k]] <- matrix(0, nrow=P+1, ncol=m)
  }
  
  ## Sylvester's equation ##
  A_array <- array(rep(0, (P+1)*(P+1)*K), c(P+1, P+1, K))
  A <- lapply(seq(dim(A_array)[3]), function(x) A_array[ , , x])
  for(k in 1:K){
    for(i in 1:n){
      A[[k]] <- A[[k]] + (1/n)*w[[l+1]][[k]][i]*(X[i,]%*%t(X[i,]))
    }
  }
  
  E <- list()
  for(k in 1:K){
    E[[k]] <- -rho * sigma[[l]][[k]]
  }
  
  W_array <- array(rep(0, (P+1)*m*K), c(P+1, m, K))
  W <- lapply(seq(dim(W_array)[3]), function(x) W_array[ , , x]) 
  for(k in 1:K){
    for(i in 1:n){
      W[[k]] <- W[[k]] + (1/n)*w[[l+1]][[k]][i]*(X[i,]%*%t(Y[i,]))
    }
  }
  
  G <- list()
  for(k in 1:K){
    G[[k]] <- W[[k]] + rho*(Ck[[1]][[k]]-U1[[1]][[k]])%*%sigma[[l]][[k]]
  }
  
  delta <- list()
  for(k in 1:K){
    delta[[k]] <- charpoly(A[[k]])
  }
  
  fAE_array <- array(rep(0, m*m*K), c(m, m, K))
  fAE <- lapply(seq(dim(fAE_array)[3]), function(x) fAE_array[ , , x])
  for(k in 1:K){
    for(j in 1:(P+2)){
      fAE[[k]] <- fAE[[k]] + delta[[k]][j]*(E[[k]]%^%(P+2-j))
    }
  }
  
  sumAGE <- as.list(array(rep(0, K), c(1, 1, K)))
  for(k in 1:K){
    for(p in 1:(P+1)){
      for(q in 0:(p-1)){
        if(p-q-1 >= 0){
          sumAGE[[k]] <- sumAGE[[k]] + delta[[k]][P+2-p]*((A[[k]]%^%(q))%*%G[[k]]%*%(E[[k]]%^%(p-q-1)))
        }
      }
    }
  }
  
  ## Update Bk ##
  Bk_sol[[2]] <- list()
  for(k in 1:K){
    Bk_sol[[2]][[k]] <- -sumAGE[[k]]%*%ginv(fAE[[k]])
  }
  
  ## Update Ck ##
  Ck[[2]] <- list()
  Hk[[2]] <- list()
  for(k in 1:K){
    Hk[[2]][[k]] <- alpha*Bk_sol[[2]][[k]] + (1-alpha)*Ck[[1]][[k]]
  }
  for(k in 1:K){
    Ck[[2]][[k]] <- MCP_thresholding(Hk[[2]][[k]]+U1[[1]][[k]], pi_[[l+1]][k], lambda[k], rho, a)
  }
  
  ## Update U1 ##
  U1[[2]] <- list()
  for(k in 1:K){
    U1[[2]][[k]] <- U1[[1]][[k]] + Hk[[2]][[k]] - Ck[[2]][[k]]
  } 
  
  ## Residuals ##
  R[[2]] <- list()
  S[[2]] <- list()
  for(k in 1:K){
    R[[2]][[k]] <- primal_residual(Bk_sol[[2]][[k]], Ck[[2]][[k]])
    S[[2]][[k]] <- dual_residual(U1[[2]][[k]], U1[[1]][[k]], rho)
  }
  
  epsilon[[2]] <- matrix(nrow=K, ncol=2)
  for(k in 1:K){
    e_pri <- epsilon_primal(Bk_sol[[2]][[k]], Ck[[2]][[k]], e_abs,e_rel)
    e_dual <- epsilon_dual(U1[[2]][[k]], e_abs, e_rel)
    epsilon[[2]][k,] <- c(e_pri, e_dual)
  }
  
  ## The iterations of updating Bk, Ck, U1 ##
  t <- 2
  while(stop_admm_lasso(R[[t]], S[[t]], epsilon[[t]])){
    
    for(k in 1:K){
      G[[k]] <- W[[k]] + rho*(Ck[[t]][[k]]-U1[[t]][[k]])%*%sigma[[l]][[k]]
    }
    
    sumAGE <- as.list(array(rep(0, K), c(1, 1, K)))
    for(k in 1:K){
      for(p in 1:(P+1)){
        for(q in 0:(p-1)){
          if(p-q-1 >= 0){
            sumAGE[[k]] <- sumAGE[[k]] + 
              delta[[k]][P+2-p]*((A[[k]]%^%(q))%*%G[[k]]%*%(E[[k]]%^%(p-q-1)))
          }
        }
      }
    }
    
    ## Update Bk ##
    Bk_sol[[t+1]] <- list()
    for(k in 1:K){
      Bk_sol[[t+1]][[k]] <- -sumAGE[[k]]%*%ginv(fAE[[k]])
    }
    
    ## Update Ck ##
    Ck[[t+1]] <- list()
    Hk[[t+1]] <- list()
    for(k in 1:K){
      Hk[[t+1]][[k]] <- alpha*Bk_sol[[t+1]][[k]] + (1-alpha)*Ck[[t]][[k]]
    }
    for(k in 1:K){
      Ck[[t+1]][[k]] <- MCP_thresholding(Hk[[t+1]][[k]]+U1[[t]][[k]], pi_[[l+1]][k], lambda[k], rho, a)
    }
    
    ## Update U1 ##
    U1[[t+1]] <- list()
    for(k in 1:K){
      U1[[t+1]][[k]] <- U1[[t]][[k]] + Hk[[t+1]][[k]] - Ck[[t+1]][[k]]
    } 
    
    ## Residuals ##
    R[[t+1]] <- list()
    S[[t+1]] <- list()
    for(k in 1:K){
      R[[t+1]][[k]] <- primal_residual(Bk_sol[[t+1]][[k]], Ck[[t+1]][[k]])
      S[[t+1]][[k]] <- dual_residual(U1[[t+1]][[k]], U1[[t]][[k]], rho)
    }
    
    epsilon[[t+1]] <- matrix(nrow=K, ncol=2)
    for(k in 1:K){
      e_pri <- epsilon_primal(Bk_sol[[t+1]][[k]], Ck[[t+1]][[k]], e_abs, e_rel)
      e_dual <- epsilon_dual(U1[[t+1]][[k]], e_abs, e_rel)
      epsilon[[t+1]][k,] <- c(e_pri, e_dual)
    }
    
    if(t == 10000){
      break
    }else{
      t <- t+1
    }
  }
  return(Bk_sol[[t]])
}


### ADMM stopping criteria ####
primal_residual <- function(B, D) B-D

dual_residual <- function(U, u, r) -r*(U-u)

epsilon_primal <- function(B, C, epsilon_absolute, epsilon_relative){
  e_abs_term <- sqrt(length(B))*epsilon_absolute
  e_rel_term <- max(norm(B,"F"), norm(C,"F"))*epsilon_relative
  return(e_abs_term + e_rel_term)
}

epsilon_dual <- function(U, epsilon_absolute, epsilon_relative){
  e_abs_term <- sqrt(length(U))*epsilon_absolute
  e_rel_term <- norm(U,"F")*epsilon_relative
  return(e_abs_term + e_rel_term)
}

stop_admm_lasso <- function(R,S,E){
  R_norm <- sapply(R, function(x) norm(x,"F"))
  S_norm <- sapply(S, function(x) norm(x,"F"))
  return(any(any(R_norm > E[,1]), any(S_norm > E[,2]))) # E[,1] is a vector of primal residuals and E[,2] is a vector of dual residuals.
}


### Thresholding functions ###
LASSO_thresholding <- function(z,p,lam,rh){
  Ck_new <- matrix(nrow=(P+1), ncol=m)
  for(r in 1:(P+1)){
    for(s in 1:m){
      if(abs(z[r,s]) <= p*lam/rh){
        Ck_new[r,s] <- 0
      }else{
        Ck_new[r,s] <- z[r,s] - sign(z[r,s])*p*lam/rh
      }
    }
  }
  return(Ck_new)
}

SCAD_thresholding <- function(z,p,lam,rh,a){
  Ck_new <- matrix(nrow=(P+1), ncol=m)
  for(r in 1:(P+1)){
    for(s in 1:m){
      if(abs(z[r,s]) <= (1+p/rh)*lam){
        if(abs(z[r,s]) <= p*lam/rh){
          Ck_new[r,s] <- 0
        }else{
          Ck_new[r,s] <- z[r,s] - sign(z[r,s])*p*lam/rh
        }
      }else if(abs(z[r,s]) >= a*lam){
        Ck_new[r,s] <- z[r,s]
      }else{
        Ck_new[r,s] <- (z[r,s] - sign(z[r,s])*a*p*lam/(rh*(a-1))) / (1-p/(rh*(a-1)))
      }
    }
  }
  return(Ck_new)
}

MCP_thresholding <- function(z,p,lam,rh,a){
  Ck_new <- matrix(nrow=(P+1), ncol=m)
  for(r in 1:(P+1)){
    for(s in 1:m){
      if(abs(z[r,s]) <= a*lam){
        if(abs(z[r,s]) <= p*lam/rh){
          Ck_new[r,s] <- 0
        }else{
          Ck_new[r,s] <- (z[r,s] - sign(z[r,s])*p*lam/rh) / (1-p/(rh*a))
        }
      }else{
        Ck_new[r,s] <- z[r,s]
      }
    }
  }
  return(Ck_new)
}


### EM algorithm ###
## E step ##
e_step_w <- function(l){
  w[[l+1]] <- list()
  wk <- matrix(nrow=n, ncol=K)
  for(k in 1:K){
    wk[,k] <- pi_[[l]][k] * density[[l]][[k]]
  }
  for(k in 1:K){
    w[[l+1]][[k]] <- prop.table(wk, margin=1)[,k]
  }
  return(w[[l+1]])
}

## Update pi_k in M step ##
m_step_pi <- function(l){
  pi_[[l+1]] <- vector()
  sum_w <- sapply(w[[l+1]], sum)
  pi_[[l+1]] <- (1/n) * sum_w
  return(pi_[[l+1]])
}

## Update B_k in M step ##
m_step_Bk_mvFMR <- function(l){
  Bk[[l+1]] <- list()
  
  WXX_array <- array(rep(0, (P+1)*(P+1)*K), c(P+1, P+1, K))
  WXX <- lapply(seq(dim(WXX_array)[3]), function(x) WXX_array[ , , x])
  for(k in 1:K){
    for(i in 1:n){
      WXX[[k]] <- WXX[[k]] + w[[l+1]][[k]][i]*(X[i,]%*%t(X[i,]))
    }
  }
  
  WXY_array <- array(rep(0, (P+1)*m*K), c(P+1, m, K))
  WXY <- lapply(seq(dim(WXY_array)[3]), function(x) WXY_array[ , , x])
  for(k in 1:K){
    for(i in 1:n){
      WXY[[k]] <- WXY[[k]] + w[[l+1]][[k]][i]*(X[i,]%*%t(Y[i,]))
    }
  }
  
  for(k in 1:K){
    Bk[[l+1]][[k]] <- ginv(WXX[[k]])%*%WXY[[k]]
  }
  
  return(Bk[[l+1]])
}

m_step_Bk_mvFMR_L <- function(l){
  Bk[[l+1]] <- list()
  Bk[[l+1]] <- admm_Bk_LASSO(l)
  
  return(Bk[[l+1]])
}

m_step_Bk_mvFMR_S <- function(l){
  Bk[[l+1]] <- list()
  Bk[[l+1]] <- admm_Bk_SCAD(l)
  
  return(Bk[[l+1]])
}

m_step_Bk_mvFMR_M <- function(l){
  Bk[[l+1]] <- list()
  Bk[[l+1]] <- admm_Bk_MCP(l)
  
  return(Bk[[l+1]])
}

m_step_Bk_oracle <- function(l){
  Bk[[l+1]] <- list()
  
  WXX_array <- array(rep(0, (P+1)*(P+1)*K), c(P+1, P+1, K))
  WXX <- lapply(seq(dim(WXX_array)[3]), function(x) WXX_array[ , , x])
  for(k in 1:K){
    for(i in 1:n){
      WXX[[k]] <- WXX[[k]] + w[[l+1]][[k]][i]*(X[i,]%*%t(X[i,]))
    }
  }
  
  WXY_array <- array(rep(0, (P+1)*m*K), c(P+1, m, K))
  WXY <- lapply(seq(dim(WXY_array)[3]), function(x) WXY_array[ , , x])
  for(k in 1:K){
    for(i in 1:n){
      WXY[[k]] <- WXY[[k]] + w[[l+1]][[k]][i]*(X[i,]%*%t(Y[i,]))
    }
  }
  
  Bk_oracle <- list()
  for(k in 1:K){
    Bk[[l+1]][[k]] <- ginv(WXX[[k]])%*%WXY[[k]]
    Bk_oracle[[k]] <- Bk[[l+1]][[k]]
    Bk_oracle[[k]][which(true_Bk[[k]] == 0)] <- 0
    Bk[[l+1]][[k]] <- Bk_oracle[[k]]
  }
  
  return(Bk[[l+1]])
}

## Update Sigma_k in M step ##
m_step_sigma <- function(l){
  sigma[[l+1]] <- list()
  sigma_sum_array <- array(rep(0, m*m*K), c(m, m, K))
  sigma_sum <- lapply(seq(dim(sigma_sum_array)[3]), function(x) sigma_sum_array[ , , x])
  
  for(k in 1:K){
    for(i in 1:n){
      sigma_sum[[k]] <- sigma_sum[[k]] + w[[l+1]][[k]][i] * 
        ((Y[i,]-t(Bk[[l+1]][[k]])%*%X[i,]) %*% t(Y[i,]-t(Bk[[l+1]][[k]])%*%X[i,]))
    }
    sigma[[l+1]][[k]] <- (1/sum(w[[l+1]][[k]])) * sigma_sum[[k]]
  }
  return(sigma[[l+1]])
}

## MSE ##
total_diff_norm <- function(l){
  pi_norm <- vector()
  Bk_norm <- vector()
  sigma_norm <- vector()
  
  for(k in 1:K){
    pi_norm[k] <- (pi_[[l+1]][k]-pi_[[l]][k])^2
    Bk_norm[k] <- sum((Bk[[l+1]][[k]]-Bk[[l]][[k]])^2) # Frobenius norm
    sigma_norm[k] <- sum((sigma[[l+1]][[k]]-sigma[[l]][[k]])^2) # Frobenius norm
  }
  return(sqrt(sum(pi_norm,Bk_norm,sigma_norm)))
}

## Round down to five decimal places ##
floor_5 <- function(x) as.data.frame(trunc(x*10^4)/10^4)

## Round down to six decimal places ##
floor_6 <- function(x) as.data.frame(trunc(x*10^5)/10^5)

## The modified BIC rounding down to five decimal places ##
BIC_5 <- function(l){
  logpi_term <- sum(sapply(w[[l]], sum) * log(pi_[[l]]))
  logdensity_term <- 0
  for(k in 1:K){
    logdensity_term <- logdensity_term + sum(w[[l]][[k]] * lapply(density[[l]], log)[[k]])
  }
  negative_ll <- (-2) * (logpi_term + logdensity_term)
  d_e <- K + (K-1) + sum(sapply(lapply(Bk[[l]], floor_5), function(x) colSums(x!=0)))
  
  return(negative_ll + log(n)*d_e)
}

## The modified BIC rounding down to six decimal places ##
BIC_6 <- function(l){
  logpi_term <- sum(sapply(w[[l]], sum) * log(pi_[[l]]))
  logdensity_term <- 0
  for(k in 1:K){
    logdensity_term <- logdensity_term + sum(w[[l]][[k]] * lapply(density[[l]], log)[[k]])
  }
  negative_ll <- (-2) * (logpi_term + logdensity_term)
  d_e <- K + (K-1) + sum(sapply(lapply(Bk[[l]], floor_6), function(x) colSums(x!=0)))
  
  return(negative_ll + log(n)*d_e)
}

