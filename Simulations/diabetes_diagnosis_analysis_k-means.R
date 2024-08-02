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
# by using k-means clustering to generate Figure 1.

# The R-code compares the results of clustering using mvFMR-SCAD and 
# k-means clustering method, when K is 3.


library(ggplot2)
library(patchwork)

source("./Data/diabetes_diagnosis_data.R")

set.seed(1032)
group <- kmeans(Y, centers=3)$cluster

dat_Y <- as.data.frame(cbind(Y, group=group))
for(i in 1:nrow(Y)){
  if(dat_Y$group[i] == 1){
    dat_Y$group[i] <- "Group 1"
  }else if(dat_Y$group[i] == 2){
    dat_Y$group[i] <- "Group 3"
  }else{
    dat_Y$group[i] <- "Group 2"
  }
}

k_means_plot <- ggplot(data=dat_Y, aes(x=glucose, y=HbA1c, group=group)) +
  geom_point(aes(shape=group, color=group, size=group)) + 
  ggthemes::theme_few() +
  theme(legend.title=element_blank())  +
  labs(x="glucose", y="HbA1c", title="k-means clustering") +
  scale_color_manual(values=c("tomato2","deepskyblue2","black","green3")) +
  scale_size_manual(values=c(2,2,2,2)) +
  scale_shape_manual(values=c(16,17,15,8)) +
  theme(axis.title.x=element_text(size=13),
        axis.title.y=element_text(size=13)) + 
  theme(plot.title=element_text(hjust=0.5)) +
  theme(legend.text=element_text(size=13)) +
  scale_x_continuous(breaks=c(100,126,250,400)) +
  scale_y_continuous(breaks=c(4,5.7,6.5,10,16))

k_means_plot2 <- k_means_plot + geom_hline(yintercept=6.5, linetype="dashed", size=0.6, col="grey50") + 
  geom_hline(yintercept=5.7,  linetype="dashed", size=0.6, col="grey50") +
  geom_vline(xintercept=126, linetype="dashed", size=0.6, col="grey50") +
  geom_vline(xintercept=100, linetype="dashed", size=0.6, col="grey50")

k_means_plot2

## Showing plots side by side ##
mvFMR_S_plot2 + k_means_plot2
