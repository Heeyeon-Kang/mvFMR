##############################################################################
####              Penalized estimation for a finite mixture of            ####
####                   multivariate regression models                     ####
##############################################################################
####                        Author: Heeyeon Kang                          ####
####                      Supervisor: Sunyoung Shin                       ####
##############################################################################
####       R-code of the data of simulations presented in Section 5.      ####
##############################################################################

# The following R-code generates the data of simulations presented in 
# Section 5.1, Section 5.2, and Section 5.3.

# We set the seed numbers to fix all the dataset used in our simulations.
# The seed numbers of Section 5.1, Section 5.2, and Section 5.3 are 
# 2,5,10 / 210,220,230 / 302,304, respectively, for each model.

# The data size is written in parentheses in the comments. 
# For generating 1000 or 2000 samples, please replace the nrow(X_data) with 1000 or 2000.


### Section 5.1 : model 1 (500) ###
###       m=2, P=10, K=2        ###
data_generate_5.1.1 <- function(z){
  
  ## The mixing proportions ##
  true_pi_ <- vector()
  true_pi_[1] <- 0.5
  true_pi_[2] <- 1 - true_pi_[1]

  ## The variance-covariance matrices of errors ##
  true_sigma <- list()
  true_sigma[[1]] <- matrix(rep(0.5, 2^2), nrow=2) + 0.5*diag(2)
  true_sigma[[2]] <- matrix(rep(0.5, 2^2), nrow=2) + 0.5*diag(2)
  
  ## The coefficients matrices of each components ##
  true_Bk <- list()
  true_Bk[[1]] <- matrix(c( 1,-1,
                            0, 2,
                            0, 0,
                            3, 0,
                           -2, 3,
                            0, 0,
                            0, 0,
                            0, 0,
                            0, 0,
                            0, 0), ncol=2, byrow=TRUE)
  true_Bk[[2]] <- matrix(c(-1, 1,
                            3, 0,
                            0, 2,
                            1,-1,
                            0, 0,
                            0, 0,
                            0, 0,
                            0, 0,
                            0, 0,
                            0, 0), ncol=2, byrow=TRUE)
  true <- list(pi=true_pi_, Bk=true_Bk, sigma=true_sigma)
  
  set.seed(z)
  
  X_data <- matrix(0, nrow=500, ncol=nrow(true_Bk[[1]])) 
  for(i in 1:nrow(X_data)){
    x_data <- rnorm(n=(ncol(X_data)-1), mean=0, sd=1)
    X_data[i,] <- cbind(1,t(x_data))
  }
  
  error <- list()
  error[[1]] <- MASS::mvrnorm(n=nrow(X_data), rep(0,nrow(true_sigma[[1]])), true_sigma[[1]])
  error[[2]] <- MASS::mvrnorm(n=nrow(X_data), rep(0,nrow(true_sigma[[1]])), true_sigma[[2]])
  
  cluster <- sample(c(1,2), size=nrow(X_data), replace=TRUE, prob=true_pi_)
  
  rb1 <- rep(0, nrow(X_data))
  rb1[which(cluster == 1)] <- 1
  rb2 <- rep(0, nrow(X_data))
  rb2[which(cluster == 2)] <- 1
  
  Y_data <- rb1*(X_data%*%true_Bk[[1]]+error[[1]]) + rb2*(X_data%*%true_Bk[[2]]+error[[2]])
  
  return(list(X=X_data, Y=Y_data, cluster=cluster, true=true))
}

### Section 5.1 : model 2 (500) ###
###       m=5, P=10, K=2        ###
data_generate_5.1.2 <- function(z){
  
  ## The mixing proportions ##
  true_pi_ <- vector()
  true_pi_[1] <- 0.5
  true_pi_[2] <- 1 - true_pi_[1]
  
  ## The variance-covariance matrices of errors ##
  true_sigma <- list()
  true_sigma[[1]] <- matrix(rep(0.5, 5^2), nrow=5) + 0.5*diag(5)
  true_sigma[[2]] <- matrix(rep(0.5, 5^2), nrow=5) + 0.5*diag(5)
  
  ## The coefficients matrices of each components ##
  true_Bk <- list()
  true_Bk[[1]] <- matrix(c( 1,-1, 1, 2,-1,
                            0, 2, 0, 1, 0,
                            0, 0,-2,-1, 3,
                            3, 0, 0, 0, 1,
                           -2, 3, 0, 1, 0,
                            0, 0, 0, 0, 0,
                            0, 0, 0, 0, 0,
                            0, 0, 0, 0, 0,
                            0, 0, 0, 0, 0,
                            0, 0, 0, 0, 0 ), ncol=5, byrow=TRUE)
  true_Bk[[2]] <- matrix(c(-1, 1, 2,-2, 1,
                            3, 0, 0, 1, 0,
                            0, 2, 1, 0, 1,
                            1,-1, 0, 0, 0,
                            0, 0, 3, 0, 0,
                            0, 0, 0, 0, 0,
                            0, 0, 0, 0, 0,
                            0, 0, 0, 0, 0,
                            0, 0, 0, 0, 0,
                            0, 0, 0, 0, 0 ), ncol=5, byrow=TRUE)
  true <- list(pi=true_pi_, Bk=true_Bk, sigma=true_sigma)
  
  set.seed(z)
  
  X_data <- matrix(0, nrow=500, ncol=nrow(true_Bk[[1]]))
  for(i in 1:nrow(X_data)){
    x_data <- rnorm(n=(ncol(X_data)-1), mean=0, sd=1)
    X_data[i,] <- cbind(1,t(x_data))
  }
  
  error <- list()
  error[[1]] <- MASS::mvrnorm(n = nrow(X_data), rep(0, nrow(true_sigma[[1]])), true_sigma[[1]])
  error[[2]] <- MASS::mvrnorm(n = nrow(X_data), rep(0, nrow(true_sigma[[1]])), true_sigma[[2]])
  
  cluster <- sample(c(1,2), size=nrow(X_data), replace=TRUE, prob=true_pi_)
  
  rb1 <- rep(0, nrow(X_data))
  rb1[which(cluster == 1)] <- 1
  rb2 <- rep(0, nrow(X_data))
  rb2[which(cluster == 2)] <- 1
  
  Y_data <- rb1*(X_data%*%true_Bk[[1]]+error[[1]]) + rb2*(X_data%*%true_Bk[[2]]+error[[2]])
  
  return(list(X=X_data, Y=Y_data, cluster=cluster, true=true))
}

### Section 5.1 : model 3 (1000) ###
###       m=10, P=10, K=2        ###
data_generate_5.1.3 <- function(z){
  
  ## The mixing proportions ##
  true_pi_ <- vector()
  true_pi_[1] <- 0.5
  true_pi_[2] <- 1 - true_pi_[1]
  
  ## The variance-covariance matrices of errors ##
  true_sigma <- list(matrix(), matrix())
  true_sigma[[1]] <- matrix(rep(0.5, 10^2), nrow=10) + 0.5*diag(10)
  true_sigma[[2]] <- matrix(rep(0.5, 10^2), nrow=10) + 0.5*diag(10)
  
  ## The coefficients matrices of each components ##
  true_Bk <- list(matrix(), matrix())
  true_Bk[[1]] <- matrix(c( 1,-1, 1, 2,-1,-2, 1, 2,-3,-1,
                            0, 2, 0, 1, 0, 0,-3, 0, 1, 0,
                            0, 0,-2,-1, 3, 1,-1, 0,-1, 0,
                            3, 0, 0, 0, 1,-1, 2, 0, 0, 1,
                           -2, 3, 0, 1, 0, 0, 0,-1, 2,-3,
                            0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                            0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                            0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                            0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                            0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ), ncol=10, byrow=TRUE)
  true_Bk[[2]] <- matrix(c(-1, 1, 2,-2, 1, 3,-3, 2,-1,-2,
                            3, 0, 0, 1, 0, 0, 2, 0, 2, 0,
                            0, 2, 1, 0, 1, 0,-1, 0, 3, 0,
                            1,-1, 0, 0, 0,-1, 1,-1, 0, 0,
                            0, 0, 3, 0, 0, 0, 0, 3, 0, 1,
                            0, 0, 0, 0, 0, 0, 0, 0, 0,-1,
                            0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                            0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                            0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                            0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ), ncol=10, byrow=TRUE)
  true <- list(pi=true_pi_, Bk=true_Bk, sigma=true_sigma)
  
  set.seed(z)
  
  X_data <- matrix(0, nrow=1000, ncol=nrow(true_Bk[[1]]))
  for(i in 1:nrow(X_data)){
    x_data <- rnorm(n=(ncol(X_data)-1), mean=0, sd=1)
    X_data[i,] <- cbind(1,t(x_data))
  }
  
  error <- list(matrix(0,nrow(X_data),ncol(true_Bk[[1]])), matrix(0,nrow(X_data),ncol(true_Bk[[1]])))
  error[[1]] <- MASS::mvrnorm(n = nrow(X_data), rep(0, nrow(true_sigma[[1]])), true_sigma[[1]])
  error[[2]] <- MASS::mvrnorm(n = nrow(X_data), rep(0, nrow(true_sigma[[1]])), true_sigma[[2]])
  
  cluster <- sample(c(1, 2), size=nrow(X_data), replace=TRUE, prob=c(true_pi_[1], true_pi_[2]))
  
  rb1 <- rep(0, nrow(X_data))
  rb1[which(cluster == 1)] <- 1
  rb2 <- rep(0, nrow(X_data))
  rb2[which(cluster == 2)] <- 1
  
  Y_data <- rb1*(X_data%*%true_Bk[[1]]+error[[1]]) + rb2*(X_data%*%true_Bk[[2]]+error[[2]])
  
  return(list(X=X_data, Y=Y_data, cluster=cluster, true=true))
}

### Section 5.2 : model 1 (500) ###
###       m=2, P=10, K=2        ###
data_generate_5.2.1 <- function(z){
  
  ## The mixing proportions ##
  true_pi_ <- vector()
  true_pi_[1] <- 0.5
  true_pi_[2] <- 1 - true_pi_[1]
  
  ## The variance-covariance matrices of errors ##
  true_sigma <- list()
  true_sigma[[1]] <- matrix(rep(0.5, 2^2), nrow=2) + 0.5*diag(2)
  true_sigma[[2]] <- matrix(rep(0.5, 2^2), nrow=2) + 0.5*diag(2)
  
  ## The coefficients matrices of each components ##
  true_Bk <- list()
  true_Bk[[1]] <- matrix(c(   1,-0.5,
                           -0.5,   0,
                              0, 1.5,
                            1.5, 0.5,
                              0,   0,
                              0,   0,
                              0,   0,
                              0,   0,
                           -0.5,   0,
                              0,-1.5), ncol=2, byrow=TRUE)
  true_Bk[[2]] <- matrix(c(  -1, 0.5,
                              0, 1.5,
                            0.5,-0.5,
                           -1.5,   0,
                              0,   1,
                              0,   0,
                              0,   0,
                              0,   0,
                            1.5,   0,
                           -0.5,   1), ncol=2, byrow=TRUE)
  true <- list(pi=true_pi_, Bk=true_Bk, sigma=true_sigma)

  set.seed(z)
  
  X_data <- matrix(0, nrow=500, ncol=nrow(true_Bk[[1]]))
  for(i in 1:nrow(X_data)){
    x_data <- rnorm(n=(ncol(X_data)-1), mean=0, sd=1)
    X_data[i,] <- cbind(1,t(x_data))
  }
  
  error <- list()
  error[[1]] <- MASS::mvrnorm(n=nrow(X_data), rep(0,nrow(true_sigma[[1]])), true_sigma[[1]])
  error[[2]] <- MASS::mvrnorm(n=nrow(X_data), rep(0,nrow(true_sigma[[1]])), true_sigma[[2]])
  
  cluster <- sample(c(1,2), size=nrow(X_data), replace=TRUE, prob=true_pi_)
  
  rb1 <- rep(0, nrow(X_data))
  rb1[which(cluster == 1)] <- 1
  rb2 <- rep(0, nrow(X_data))
  rb2[which(cluster == 2)] <- 1
  
  Y_data <- rb1*(X_data%*%true_Bk[[1]]+error[[1]]) + rb2*(X_data%*%true_Bk[[2]]+error[[2]])
  
  return(list(X=X_data, Y=Y_data, cluster=cluster, true=true))
}

### Section 5.2 : model 2 (500) ###
###       m=2, P=20, K=2        ###
data_generate_5.2.2 <- function(z){
  
  ## The mixing proportions ##
  true_pi_ <- vector()
  true_pi_[1] <- 0.5
  true_pi_[2] <- 1 - true_pi_[1]
  
  ## The variance-covariance matrices of errors ##
  true_sigma <- list()
  true_sigma[[1]] <- matrix(rep(0.5, 2^2), nrow=2) + 0.5*diag(2)
  true_sigma[[2]] <- matrix(rep(0.5, 2^2), nrow=2) + 0.5*diag(2)
  
  ## The coefficients matrices of each components ##
  true_Bk <- list()
  zero_mat <- matrix(0, nrow=10, ncol=2)
  true_Bk[[1]] <- matrix(c(   1,-0.5,
                           -0.5,   0,
                              0, 1.5,
                            1.5, 0.5,
                              0,   0,
                              0,   0,
                              0,   0,
                              0,   0,
                           -0.5,   0,
                              0,-1.5), ncol=2, byrow=TRUE)
  true_Bk[[2]] <- matrix(c(  -1, 0.5,
                              0, 1.5,
                            0.5,-0.5,
                           -1.5,   0,
                              0,   1,
                              0,   0,
                              0,   0,
                              0,   0,
                            1.5,   0,
                           -0.5,   1), ncol=2, byrow=TRUE)
  true_Bk <- lapply(true_Bk, function(x) rbind(x, zero_mat))
  true <- list(pi=true_pi_, Bk=true_Bk, sigma=true_sigma)
  
  set.seed(z)
  
  X_data <- matrix(0, nrow=500, ncol=nrow(true_Bk[[1]]))
  for(i in 1:nrow(X_data)){
    x_data <- rnorm(n=(ncol(X_data)-1), mean=0, sd=1)
    X_data[i,] <- cbind(1,t(x_data))
  }
  
  error <- list()
  error[[1]] <- MASS::mvrnorm(n=nrow(X_data), rep(0,nrow(true_sigma[[1]])), true_sigma[[1]])
  error[[2]] <- MASS::mvrnorm(n=nrow(X_data), rep(0,nrow(true_sigma[[2]])), true_sigma[[2]])
  
  cluster <- sample(c(1,2), size=nrow(X_data), replace=TRUE, prob=true_pi_)
  
  rb1 <- rep(0, nrow(X_data))
  rb1[which(cluster == 1)] <- 1
  rb2 <- rep(0, nrow(X_data))
  rb2[which(cluster == 2)] <- 1
  
  Y_data <- rb1*(X_data%*%true_Bk[[1]]+error[[1]]) + rb2*(X_data%*%true_Bk[[2]]+error[[2]])
  
  return(list(X=X_data, Y=Y_data, cluster=cluster, true=true))
}

### Section 5.2 : model 3 (500) ###
###       m=2, P=30, K=2        ###
data_generate_5.2.3 <- function(z){
  
  ## The mixing proportions ##
  true_pi_ <- vector()
  true_pi_[1] <- 0.5
  true_pi_[2] <- 1 - true_pi_[1]
  
  ## The variance-covariance matrices of errors ##
  true_sigma <- list(matrix(), matrix())
  true_sigma[[1]] <- matrix(rep(0.5, 2^2), nrow=2) + 0.5*diag(2)
  true_sigma[[2]] <- matrix(rep(0.5, 2^2), nrow=2) + 0.5*diag(2)
  
  ## The coefficients matrices of each components ##
  true_Bk <- list(matrix(), matrix())
  zero_mat <- matrix(0, nrow=20, ncol=2)
  true_Bk[[1]] <- matrix(c(   1,-0.5,
                           -0.5,   0,
                              0, 1.5,
                            1.5, 0.5,
                              0,   0,
                              0,   0,
                              0,   0,
                              0,   0,
                           -0.5,   0,
                              0,-1.5), ncol=2, byrow=TRUE)
  true_Bk[[2]] <- matrix(c(  -1, 0.5,
                              0, 1.5,
                            0.5,-0.5,
                           -1.5,   0,
                              0,   1,
                              0,   0,
                              0,   0,
                              0,   0,
                            1.5,   0,
                           -0.5,   1), ncol=2, byrow=TRUE)
  true_Bk <- lapply(true_Bk, function(x) rbind(x, zero_mat))
  true <- list(pi=true_pi_, Bk=true_Bk, sigma=true_sigma)
  
  set.seed(z)
  
  X_data <- matrix(0, nrow=500, ncol=nrow(true_Bk[[1]]))
  for(i in 1:nrow(X_data)){
    x_data <- rnorm(n=(ncol(X_data)-1), mean=0, sd=1)
    X_data[i,] <- cbind(1,t(x_data))
  }
  
  error <- list()
  error[[1]] <- MASS::mvrnorm(n=nrow(X_data), rep(0,nrow(true_sigma[[1]])), true_sigma[[1]])
  error[[2]] <- MASS::mvrnorm(n=nrow(X_data), rep(0,nrow(true_sigma[[2]])), true_sigma[[2]])
  
  cluster <- sample(c(1,2), size=nrow(X_data), replace=TRUE, prob=true_pi_)
  
  rb1 <- rep(0, nrow(X_data))
  rb1[which(cluster == 1)] <- 1
  rb2 <- rep(0, nrow(X_data))
  rb2[which(cluster == 2)] <- 1
  
  Y_data <- rb1*(X_data%*%true_Bk[[1]]+error[[1]]) + rb2*(X_data%*%true_Bk[[2]]+error[[2]])
  
  return(list(X=X_data, Y=Y_data, cluster=cluster, true=true))
}

### Section 5.3 : model 1 (500) ###
###       m=2, P=10, K=2        ###
data_generate_5.3.1 <- function(z){
  
  ## The mixing proportions ##
  true_pi_ <- vector()
  true_pi_[1] <- 0.5
  true_pi_[2] <- 1 - true_pi_[1]
  
  ## The variance-covariance matrices of errors ##
  true_sigma <- list()
  true_sigma[[1]] <- matrix(rep(0.5, 2^2), nrow=2) + 0.5*diag(2)
  true_sigma[[2]] <- matrix(rep(0.5, 2^2), nrow=2) + 0.5*diag(2)
  
  ## The coefficients matrices of each components ##
  true_Bk <- list()
  true_Bk[[1]] <- matrix(c(-5, 2,
                            2,-4,
                            4, 3,
                            0, 0,
                            0, 0,
                            0, 0,
                            0, 0,
                            0, 0,
                            0, 0,
                            0, 0), ncol=2, byrow=TRUE)
  true_Bk[[2]] <- matrix(c( 3, 4,
                            0, 0,
                            0, 0,
                           -3, 5,
                            0, 0,
                            0, 0,
                            4,-3,
                            5, 2,
                            0, 0,
                            0, 0), ncol=2, byrow=TRUE)
  true <- list(pi=true_pi_, Bk=true_Bk, sigma=true_sigma)
  
  set.seed(z)
  
  X_data <- matrix(0, nrow=500, ncol=nrow(true_Bk[[1]])) 
  for(i in 1:nrow(X_data)){
    x_data <- rnorm(n=(ncol(X_data)-1), mean=0, sd=1)
    X_data[i,] <- cbind(1,t(x_data))
  }
  
  error <- list()
  error[[1]] <- MASS::mvrnorm(n=nrow(X_data), rep(0,nrow(true_sigma[[1]])), true_sigma[[1]])
  error[[2]] <- MASS::mvrnorm(n=nrow(X_data), rep(0,nrow(true_sigma[[2]])), true_sigma[[2]])
  
  cluster <- sample(c(1,2), size=nrow(X_data), replace=TRUE, prob=true_pi_)
  
  rb1 <- rep(0, nrow(X_data))
  rb1[which(cluster == 1)] <- 1
  rb2 <- rep(0, nrow(X_data))
  rb2[which(cluster == 2)] <- 1
  
  Y_data <- rb1*(X_data%*%true_Bk[[1]]+error[[1]]) + rb2*(X_data%*%true_Bk[[2]]+error[[2]])
  
  return(list(X=X_data, Y=Y_data, cluster=cluster, true=true))
}

### Section 5.3 : model 2 (500) ###
###       m=2, P=10, K=4        ###
data_generate_5.3.2 <- function(z){
  
  ## The mixing proportions ##
  true_pi_ <- vector()
  true_pi_[1] <- 0.25
  true_pi_[2] <- 0.25
  true_pi_[3] <- 0.25
  true_pi_[4] <- 1 - sum(true_pi_[1:3])
  
  ## The variance-covariance matrices of errors ##
  true_sigma <- list()
  true_sigma[[1]] <- matrix(rep(0.5, 2^2), nrow=2) + 0.5*diag(2)
  true_sigma[[2]] <- matrix(rep(0.5, 2^2), nrow=2) + 0.5*diag(2)
  true_sigma[[3]] <- matrix(rep(0.5, 2^2), nrow=2) + 0.5*diag(2)
  true_sigma[[4]] <- matrix(rep(0.5, 2^2), nrow=2) + 0.5*diag(2)
  
  ## The coefficients matrices of each components ##
  true_Bk <- list()
  true_Bk[[1]] <- matrix(c(-5, 2,
                            2,-4,
                            4, 3,
                            0, 0,
                            0, 0,
                            0, 0,
                            0, 0,
                            0, 0,
                            0, 0,
                            0, 0), ncol=2, byrow=TRUE)
  true_Bk[[2]] <- matrix(c( 3, 4,
                            0, 0,
                            0, 0,
                           -3, 5,
                            0, 0,
                            0, 0,
                            4,-3,
                            5, 2,
                            0, 0,
                            0, 0), ncol=2, byrow=TRUE)
  true_Bk[[3]] <- matrix(c( 2,-4,
                            0, 0,
                            0, 0,
                            0, 0,
                           -3, 3,
                            0, 0,
                            0, 0,
                            0, 0,
                            5, 2,
                            0, 0), ncol=2, byrow=TRUE)
  true_Bk[[4]] <- matrix(c(-4, 3,
                            0, 0,
                            0, 0,
                            0, 0,
                            0, 0,
                           -2, 4,
                            0, 0,
                            0, 0,
                            0, 0,
                            3,-5), ncol=2, byrow=TRUE)
  true <- list(pi=true_pi_, Bk=true_Bk, sigma=true_sigma)
  
  set.seed(z)
  
  X_data <- matrix(0, nrow=500, ncol=nrow(true_Bk[[1]])) 
  for(i in 1:nrow(X_data)){
    x_data <- rnorm(n=(ncol(X_data)-1), mean=0, sd=1)
    X_data[i,] <- cbind(1,t(x_data))
  }
  
  error <- list()
  error[[1]] <- MASS::mvrnorm(n=nrow(X_data), rep(0,nrow(true_sigma[[1]])), true_sigma[[1]])
  error[[2]] <- MASS::mvrnorm(n=nrow(X_data), rep(0,nrow(true_sigma[[2]])), true_sigma[[2]])
  error[[3]] <- MASS::mvrnorm(n=nrow(X_data), rep(0,nrow(true_sigma[[3]])), true_sigma[[3]])
  error[[4]] <- MASS::mvrnorm(n=nrow(X_data), rep(0,nrow(true_sigma[[4]])), true_sigma[[4]])
  
  cluster <- sample(c(1,2,3,4), size=nrow(X_data), replace=TRUE, prob=true_pi_)
  
  rb1 <- rep(0, nrow(X_data))
  rb1[which(cluster == 1)] <- 1
  rb2 <- rep(0, nrow(X_data))
  rb2[which(cluster == 2)] <- 1
  rb3 <- rep(0, nrow(X_data))
  rb3[which(cluster == 3)] <- 1
  rb4 <- rep(0, nrow(X_data))
  rb4[which(cluster == 4)] <- 1
  
  Y_data <- rb1*(X_data%*%true_Bk[[1]]+error[[1]]) + rb2*(X_data%*%true_Bk[[2]]+error[[2]]) + 
    rb3*(X_data%*%true_Bk[[3]]+error[[3]]) + rb4*(X_data%*%true_Bk[[4]]+error[[4]])
  
  return(list(X=X_data, Y=Y_data, cluster=cluster, true=true))
}
