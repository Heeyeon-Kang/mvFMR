##############################################################################
####              Penalized estimation for a finite mixture of            ####
####                   multivariate regression models                     ####
##############################################################################
####                        Author: Heeyeon Kang                          ####
####                      Supervisor: Sunyoung Shin                       ####
##############################################################################
####  R-code of the data of real data analysis presented in Section 6.1.  ####
##############################################################################

# The raw data can be accessed from https://hbiostat.org/data/.
# Data is cleaned by removing unnecessary variables and missing observations 
# and using dummy variables.
# Also, the units are different for each variable, 
# so we did centering and scaling.

# After data cleansing, there exists 375 samples.

#diabetes_dat <- read.csv(file="~/Desktop/diabetes_2.csv")
diabetes_dat <- read.csv(file="./data/diabetes_diagnosis.csv")

delete_columns <- c(1,7,12,15,16)
diabetes_dat <- diabetes_dat[,-delete_columns]
diabetes_dat <- na.omit(diabetes_dat)

gender <- ifelse(diabetes_dat$gender == "male", 0, 1)

diabetes_Y <- as.matrix(data.frame(glucose=diabetes_dat$stab.glu, HbA1c=diabetes_dat$glyhb))
diabetes_X_1 <- as.matrix(data.frame(chol=diabetes_dat$chol,
                                     hdl=diabetes_dat$hd,
                                     age=diabetes_dat$age))
diabetes_X_2 <- as.matrix(data.frame(height=diabetes_dat$height,
                                     weight=diabetes_dat$weight,
                                     bp.1s=diabetes_dat$bp.1s,
                                     bp.1d=diabetes_dat$bp.1d,
                                     waist=diabetes_dat$waist,
                                     hip=diabetes_dat$hip,
                                     time.ppn=diabetes_dat$time.ppn))

X_before <- cbind(intercept=rep(1, nrow(diabetes_Y)), diabetes_X_1, gender, diabetes_X_2)
X <- cbind(scale(diabetes_X_1, center=T, scale=T), gender, scale(diabetes_X_2, center=T, scale=T))
X <- cbind(intercept=rep(1, nrow(X)), X)
Y <- diabetes_Y