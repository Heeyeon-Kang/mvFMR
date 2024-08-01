##############################################################################
####              Penalized estimation for a finite mixture of            ####
####                   multivariate regression models                     ####
##############################################################################
####                        Author: Heeyeon Kang                          ####
####                      Supervisor: Sunyoung Shin                       ####
##############################################################################
####     R-code of analyzing the result of the diabetes diagnosis data    ####
####                      presented in Section 6.1.                       ####
##############################################################################

# The following R-code analyzes the results of the diabetes diagnosis data 
# and generates Table 7 and Figure 1.

# The results of the analysis of the diabetes diagnosis data are stored in 
# diabetes_diagnosis.rda, which consists of the list of output, density, w, and BIC.

# List of 4
# - output  : The modified BIC minimized by searching 39 lambdas 
#             for K = 1, 2, 3, 4, 5.
# - w       : The conditional expectation of z_{ik} at E step which is 
#             the mixing proportion of the data points and is utilized 
#             for clustering.
# - density : The multivariate Gaussian distribution corresponding to each 
#             data point.
# - BIC     : The modified BIC for K = 1, 2, 3, 4, 5.

# The value of K for which the modified BIC has the smallest value is 4.
# However, when K is 4, the mixing proportion of group 3 has a very small value of 0.006. 
# Therefore, we considered the second smallest K=3, and the modified BIC 
# is not significantly different from K=4.


library(ggplot2)

load(file="./Simulations/diabetes_diagnosis.rda")
source("./Data/diabetes_diagnosis_data.R")
source("./Functions/functions.R")

sapply(diabetes_diagnosis$BIC, min)

## Optimal K = 4 ##
optimal_K <- which.min(diabetes_diagnosis$BIC)
optimal_output <- diabetes_diagnosis$output[[optimal_K]]
optimal_w <- diabetes_diagnosis$w[[optimal_K]]
optimal_density <- diabetes_diagnosis$density[[optimal_K]]

## K = 3 ##
analysis_K <- 3
analysis_output <- diabetes_diagnosis$output[[analysis_K]]
analysis_w <- diabetes_diagnosis$w[[analysis_K]]
analysis_density <- diabetes_diagnosis$density[[analysis_K]]

### Generating the Table 7 ###
# Round down to five decimal places #
Bk_floor5 <- lapply(analysis_output$Bk, floor_5)

# To restore centering and scaling #
sd_X <- apply(X_before, 2, sd)
Bk_new <- Bk_floor5
for(k in 1:analysis_K){
  for(i in 2:ncol(X)){
    for(j in 1:ncol(Y)){
      if(Bk_floor5[[k]][i,j] != 0){
        Bk_new[[k]][i,j] <- Bk_floor5[[k]][i,j] / sd_X[i]
      }
    }
  }
}

# Calculating the correlation between glucose and HbA1c #
correlation <- sapply(analysis_output$sigma, function(x) x[1,2]/sqrt(x[1,1]*x[2,2]))

# Combining to make the Table 7 #
diabetes_diagnosis_analysis <- list(pi=analysis_output$pi, Bk=Bk_new, corr=correlation)


### Generating the Figure 1 ###
group <- vector(length=nrow(Y))
for(i in 1:length(group)){
  group[i] <- which.max(lapply(analysis_w, function(x) x[i]))
}

dat_Y <- as.data.frame(cbind(Y, group=group))
for(i in 1:length(group)){
  if(dat_Y$group[i] == 1){
    dat_Y$group[i] <- "Group 1"
  }else if(dat_Y$group[i] == 2){
    dat_Y$group[i] <- "Group 2"
  }else{
    dat_Y$group[i] <- "Group 3"
  }
}

# The scatter plot by group clustered using mvFMR_S #
mvFMR_S_plot <- ggplot(data=dat_Y, aes(x=glucose, y=HbA1c, group=group)) +
  geom_point(aes(shape=group, color=group, size=group)) + 
  ggthemes::theme_few() + 
  theme(legend.title = element_blank())  +
  labs(x = "glucose", y = "HbA1c", title = "mvFMR-SCAD") +
  scale_color_manual(values=c("tomato2","deepskyblue2","black","green3")) +
  scale_size_manual(values=c(2,2,2,2)) +
  scale_shape_manual(values=c(16,17,15,8)) +
  theme(axis.title.x = element_text(size=13),
        axis.title.y = element_text(size=13)) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.text = element_text(size=13))

# The diabetes diagnosis criteria by American Diabetes Association(ADA) #
mvFMR_S_plot2 <- mvFMR_S_plot + geom_hline(yintercept=6.5, linetype="dashed", size=0.6, col="grey50") + 
  geom_hline(yintercept=5.7, linetype="dashed", size=0.6, col="grey50") +
  geom_vline(xintercept=126, linetype="dashed", size=0.6, col="grey50") +
  geom_vline(xintercept=100, linetype="dashed", size=0.6, col="grey50")

mvFMR_S_plot2
