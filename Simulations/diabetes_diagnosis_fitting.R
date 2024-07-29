##############################################################################
####              Penalized estimation for a finite mixture of            ####
####                   multivariate regression models                     ####
##############################################################################
####                        Author: Heeyeon Kang                          ####
####                      Supervisor: Sunyoung Shin                       ####
##############################################################################
####            R-code of fitting the diabetes diagnosis data             ####
####                      presented in Section 6.1.                       ####
##############################################################################

# The following R-code fits the diabetes diagnosis data presented in Section 6.1 
# using parallel computing.

# The data is analyzed using mvFMR-SCAD and can be analyzed using the other methods
# by replacing m_step_Bk_mvFMR_S in the diabetes_diagnosis function by 
# m_step_Bk_mvFMR_L or m_step_Bk_mvFMR_M.


## Packages ##
requiredPackages <- c("MASS", "pracma", "expm", "corrcoverage", 
                      "causact", "fossil", "flexmix", "gtools",
                      "parallel", "foreach", "doParallel")
for(p in requiredPackages){
  if(!require(p,character.only = TRUE)) install.packages(p)
  library(p,character.only = TRUE)
}

## Data ##
source("./data/diabetes_diagnosis_data.R")

## Functions ##
source("./code/functions.R")

diabetes_diagnosis <- function(X, Y, K){
  
  n <- nrow(X)
  P <- ncol(X) - 1
  m <- ncol(Y)
  decimal <- 5
  
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
  
  m_step_pi <- function(l){
    pi_[[l+1]] <- vector()
    sum_w <- sapply(w[[l+1]], sum)
    pi_[[l+1]] <- (1/n) * sum_w
    return(pi_[[l+1]])
  }
  
  m_step_Bk_mvFMR_S <- function(l){
    Bk[[l+1]] <- list()
    Bk[[l+1]] <- admm_Bk_SCAD(l)
    
    return(Bk[[l+1]])
  }
  
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
  
  floor_5 <- function(x) as.data.frame(trunc(x*10^4)/10^4)
  
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
  
  lambda_s <- c(1e-4, 5*1e-3, 10*1e-3, 15*1e-3, 20*1e-3, 25*1e-3,
                30*1e-3, 35*1e-3, 40*1e-3, 45*1e-3, 50*1e-3, 55*1e-3,
                60*1e-3, 65*1e-3, 70*1e-3, 75*1e-3, 80*1e-3, 85*1e-3,
                90*1e-3, 95*1e-3, 1e-1, 15*1e-2, 2*1e-1, 25*1e-2, 3*1e-1,
                35*1e-2, 4*1e-1, 45*1e-2, 5*1e-1, 55*1e-2, 6*1e-1, 65*1e-2,
                7*1e-1, 75*1e-2, 8*1e-1, 85*1e-2, 9*1e-1, 95*1e-2, 1)
  
  
  BIC_SCAD <- vector(length=length(lambda_s))
  total_OUTPUTS_SCAD <- list(list())
  w_s <- list(list())
  density_s <- list(list())
  
  rho <- 0.1
  e_abs <- 1e-6
  e_rel <- 1e-4
  alpha <- 1.5
  a <- 10
  
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
      for(k in 1:K){
        sigma[[1]][[k]] <- k * sigma[[1]][[k]]
      }
      
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
      
      w_array <- array(rep(1/K, n*1*K), c(n, 1, K))
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
        cat("K :", K, "\n")
        cat("Iteration :", t+1, "\n")
        cat("Now, Lambda :", lambda_s[z], "\n")
        
        if(t == 100){
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
      message("Prefer 5 or 6")
    }
    
    total_OUTPUTS_SCAD[[z]] <- output
    w_s[[z]] <- w[[t]]
    density_s[[z]] <- density[[t]]
  }
  
  BIC_SCAD_f <- replace(BIC_SCAD, which(is.infinite(BIC_SCAD) == T), NA)
  optimal_OUTPUTS_SCAD <- total_OUTPUTS_SCAD[[which.min(BIC_SCAD_f)]]
  optimal_w <- w_s[[which.min(BIC_SCAD_f)]]
  optimal_density <- density_s[[which.min(BIC_SCAD_f)]]
  
  return(list(optimal=optimal_OUTPUTS_SCAD, optimal_w=optimal_w, optimal_density=optimal_density,
              total=total_OUTPUTS_SCAD, total_w=w_s, total_density=density_s,
              BIC=BIC_SCAD_f))
}

numCores <- detectCores()
registerDoParallel(numCores)

parallel_OUTPUTS <- foreach(K = 1:5,
                            .combine = list,
                            .multicombine = TRUE) %dopar% {
                              diabetes_diagnosis(X, Y, K)
                            }

BIC_SCAD_s <- lapply(parallel_OUTPUTS, function(x) x$BIC)
optimal_OUTPUTS_s <- lapply(parallel_OUTPUTS, function(x) x$optimal)
optimal_w_s <- lapply(parallel_OUTPUTS, function(x) x$optimal_w)
optimal_density_s <- lapply(parallel_OUTPUTS, function(x) x$optimal_density)
total_OUTPUTS_s <- lapply(parallel_OUTPUTS, function(x) x$total)
total_w_s <- lapply(parallel_OUTPUTS, function(x) x$total_w)
total_density_s <- lapply(parallel_OUTPUTS, function(x) x$total_density)

save(BIC_SCAD_s, file="./BIC_SCAD_diabetes_diagnosis.rdata")
save(optimal_OUTPUTS_s, file="./optimal_OUTPUTS_diabetes_diagnosis.rdata")
save(optimal_w_s, file="./optimal_w_diabetes_diagnosis.rdata")
save(optimal_density_s, file="./optimal_density_diabetes_diagnosis.rdata")
save(total_OUTPUTS_s, file="./total_OUTPUTS_diabetes_diagnosis.rdata")
save(total_w_s, file="./total_w_diabetes_diagnosis.rdata")
save(total_density_s, file="./total_density_diabetes_diagnosis.rdata")

##############################################################################
# The whole simulation was run on a Linux Cluster using R version 4.1.2.
