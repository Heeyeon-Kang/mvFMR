##############################################################################
####              Penalized estimation for a finite mixture of            ####
####                   multivariate regression models                     ####
##############################################################################
####                        Author: Heeyeon Kang                          ####
####                      Supervisor: Sunyoung Shin                       ####
##############################################################################
####  R-code of the data of real data analysis presented in Section 6.2.  ####
##############################################################################

# The raw data can be accessed from https://www.kaggle.com/datasets/harbhajansingh21/bike-sharing-dataset.
# Data is cleaned by removing unnecessary variables and missing observations 
# and using dummy variables.
# Also, the units are different for each variable, 
# so we did centering and scaling.

# After data cleansing, there exists 17397 samples.

#bike_sharing_dat <- read.csv(file="~/Desktop/hour.csv")
#bike_sharing_dat <- read.csv(file="~/Volume2/heeyeon/simulation04/model2/scripts/hour.csv")
bike_sharing_dat <- read.csv(file="./data/bike_sharing.csv")

rm_col <- c(1,2,4,5,6,8,17)
bike_dat <- bike_sharing_dat[,-rm_col]

weather_dat <- matrix(nrow=nrow(bike_dat), ncol=3)
weather_dat <- as.data.frame(weather_dat)
colnames(weather_dat) <- c("weather_2", "weather_3", "weather_4")
for(i in 1:nrow(bike_dat)){
  if(bike_dat$weathersit[i] == 1){
    weather_dat[i,] <- c(0,0,0)
  }else if(bike_dat$weathersit[i] == 2){
    weather_dat[i,] <- c(1,0,0)
  }else if(bike_dat$weathersit[i] == 3){
    weather_dat[i,] <- c(0,1,0)
  }else{
    weather_dat[i,] <- c(0,0,1)
  }
}

season_dat <- matrix(nrow=nrow(bike_dat), ncol=3)
season_dat <- as.data.frame(season_dat)
colnames(season_dat) <- c("spring", "summer", "fall")
for(i in 1:nrow(bike_dat)){
  if(bike_dat$season[i] == 1){
    season_dat[i,] <- c(0,0,0)
  }else if(bike_dat$season[i] == 2){
    season_dat[i,] <- c(1,0,0)
  }else if(bike_dat$season[i] == 3){
    season_dat[i,] <- c(0,1,0)
  }else{
    season_dat[i,] <- c(0,0,1)
  }
}

bike_dat$temp <- bike_dat$temp*41
bike_dat$atemp <- bike_dat$atemp*50
bike_dat$hum <- bike_dat$hum*100
bike_dat$windspeed <- bike_dat$windspeed*67

X_before <- cbind(intercept=rep(1, nrow(bike_dat)), season_dat, bike_dat[,2:3], weather_dat, bike_dat[,5:8])

bike_dat$temp <- scale(bike_dat$temp, center=T, scale=T)
bike_dat$atemp <- scale(bike_dat$atemp, center=T, scale=T)
bike_dat$hum <- scale(bike_dat$hum, center=T, scale=T)
bike_dat$windspeed <- scale(bike_dat$windspeed, center=T, scale=T)

X <- cbind(season_dat, bike_dat[,2:3], weather_dat, bike_dat[,5:8])
X <- as.matrix(cbind(intercept=rep(1, nrow(X)), X))
Y <- as.matrix(bike_dat[,9:10])
