#################################################
#####  Stationarization of the time series  #####
#################################################
## Prerequisites
rm(list=ls())
library(dplyr)
load("../../Dataset/bikemi_data.RData")
data <- bike_nil_week_hour
head(data)

weather_data <- read.table("../../Dataset/data_weather.csv", header = T, sep = ",")
head(weather_data)
##----------------------------------------------##

#####  Weather
## Data preprocessing
## REMARK: run only if data are not alredy preprocessed in rain
w2bin <- function(x)  ## from weather to binary
{
  return(floor(as.numeric(x)/3))
}

addon <- as.data.frame(w2bin(weather_data$Weather))
names(addon) <- "rain"

weather_data <- cbind(weather_data, addon)

write.csv(weather_data, file = "../../Dataset/data_weather.csv", sep = ",", row.names = F)
##---------------------------------------------##

## Remarkable weather plots
## Temperature
plot_temperature <- function (weather_data)
{
  plot(weather_data$Day, weather_data$Tmean, col="red", pch=19, type="l", ylim = c(min(weather_data$Tmin), max(weather_data$Tmax)+5), xlab = "Day", ylab = "Temp", main = "Temperatures")
  points(weather_data$Day, weather_data$Tmin, col="lightblue", pch=19, type="l")
  points(weather_data$Day, weather_data$Tmax, col="blue", pch=19, type="l")
  
  points(weather_data$Day, weather_data$Tmin, col="lightblue", pch=weather_data$rain+15)
  points(weather_data$Day, weather_data$Tmean, col="red", pch=weather_data$rain+15)
  points(weather_data$Day, weather_data$Tmax, col="blue", pch=weather_data$rain+15)
  
  abline(v=seq(7.5,42, by=7),col="green")
  
  legend(5, 20, legend = c("Not rainy","Rainy"), pch=c(15,16))
}

## Apparently no evident link between temperature and rain

#Humidity
plot_humidity <- function (weather_data)
{
  plot(weather_data$Day, weather_data$Humidity, col="red", type="l", xlab = "Day", ylab = "Humidity")
  points(weather_data$Day, weather_data$Humidity, col="red", pch=19)
  
  abline(v=seq(7.5,42, by=7),col="green")
}

plot_wind <- function (weather_data)
{
  plot(weather_data$Day, weather_data$Wind_mean, col="red", pch=19, type="l", ylim = c(min(weather_data$Wind_mean), max(weather_data$Wind_mean)+14), xlab = "Day", ylab = "Wind", main = "Wind")
  points(weather_data$Day, weather_data$Wind_max, col="blue", pch=19, type="l")
  
  points(weather_data$Day, weather_data$Wind_mean, col="red", pch=weather_data$rain+15)
  points(weather_data$Day, weather_data$Wind_max, col="blue", pch=weather_data$rain+15)
  
  abline(v=seq(7.5,42, by=7),col="green")
  
  legend(5, 27, legend = c("Not rainy","Rainy"), pch=c(15,16))
}

##---------------------------------------------##

#####  Daily traffic
##  Functions to test traffic
plot_day_by_day_time_series <- function(starting,ending)
{
  sub <- filter(data, StazioneInizio == starting & StazioneFine == ending)
  sub <- sub %>%
    group_by(NewDay) %>%
    summarize(n())
  names(sub) <- c("Day","N_travels")
  vec <- rep(0,42)
  vec[sub$Day] <- sub$N_travels
  
  plot(1:42,vec,pch=19,col="red", type = "l", main= paste(starting," - ",ending,sep=""), xlab= "day", ylab="counts")
  points(1:42,vec,pch=19,col="red")
  abline(v=seq(7.5,42, by=7),col="green")
}

plot_day_by_day_time_series_together <- function()
{
  sub <- data %>%
    group_by(NewDay) %>%
    summarize(n())
  names(sub) <- c("Day","N_travels")
  vec <- rep(0,42)
  vec[sub$Day] <- sub$N_travels
  
  plot(1:42,vec,pch=19,col="red", type = "l", main= "All stations together", xlab= "day", ylab="counts")
  points(1:42,vec,pch=19,col="red")
  abline(v=seq(7.5,42, by=7),col="green")
}
##---------------------------------------------##

## Testing traffic
starting <- "Cadorna 2"
ending <- "Duomo"

x11()
test0 <- plot_day_by_day_time_series(starting, ending)
x11()
test1 <- plot_day_by_day_time_series_together()

graphics.off()

##---------------------------------------------##
## Comparison traffic-weather
x11()
par(mfrow=c(1,2))
plot_day_by_day_time_series_together()
plot_temperature(weather_data)
## Apparently quite good link between temerature and bies
## Also rain is quite intersting [[might be useful to know avg mm/h]]

## Comparison traffic-humidity
x11()
par(mfrow=c(1,2))
plot_day_by_day_time_series_together()
plot_humidity(weather_data)
## Apparently irrelevant

## Comparison traffic-wind
x11()
par(mfrow=c(1,2))
plot_day_by_day_time_series_together()
plot_wind(weather_data)
# Apparently irrelevant

graphics.off()

#####---------------------------------------#####
##### Global traffic stationarization test  #####
#####---------------------------------------#####
#### Frequentist approach
d = 7 # period is a week
q = 3
##---------------------------------------------##

## Define filters
mean_filter <- function(total_traffic, q, d)
{
  tot_days <- dim(total_traffic)[1]
  mean_comp = rep(0, tot_days)
  
  for(i in 1:tot_days)
  {
    mean_comp[i] = sum(total_traffic[max(1,i-q):min(tot_days,i+q),2])/d
  }
  return(mean_comp)
}

seasonal_filter <- function(total_traffic, q, d, m)
{
  tot_days <- dim(total_traffic)[1]
  w = rep(0, d)
  
  for(k in 1:d)
  {
    w[k] = sum(total_traffic[which(as.logical(((total_traffic[,1]%%d)==k)*(k!=d)+((total_traffic[,1]%%d)==0)*(k==d))), 2])/(tot_days/d)
  }
  
  s = w - sum(w)/d
  
  periodic_comp = rep(0, tot_days)
  for (i in 1:tot_days)
  {
    periodic_comp[i] = s[(i%%d)*((i%%d)!=0)+d*((i%%d)==0)]
  }
  
  return(periodic_comp)
}
##---------------------------------------------##

## Define plotting functions
plot_mean <- function(total_traffic, m)
{
  vec <- rep(0,42)
  vec[total_traffic$Day] <- total_traffic$N_travels
  
  x11()
  par(mfrow=c(1,2))
  
  plot(1:42,vec,pch=19,col="red", type = "l", main= paste("Day vs mean"), xlab= "day", ylab="counts")
  points(1:42, m, col="blue", pch=19, type="l")
  
  points(1:42,vec,pch=19,col="red")
  points(1:42, m, col="blue", pch=19)
  
  abline(v=seq(7.5,42, by=7),col="green")
  
  plot(1:42,vec,pch=19,col="red", type = "l", main= paste("Day vs avg(mean)"), xlab= "day", ylab="counts")
  points(1:42, rep(mean(m),42), col="blue", pch=19, type="l")
  
  points(1:42,vec,pch=19,col="red")
  points(1:42, rep(mean(m),42), col="blue", pch=19)
  
  abline(v=seq(7.5,42, by=7),col="green")
}

plot_mean_periodic <- function(total_traffic, m, s)
{
  vec <- rep(0,42)
  vec[total_traffic$Day] <- total_traffic$N_travels
  
  x11()
  par(mfrow=c(1,2))
  plot(1:42,vec,pch=19,col="red", type = "l", main= paste("Day vs mean + period"), xlab= "day", ylab="counts")
  points(1:42, m+s, col="blue", pch=19, type="l")
  
  points(1:42,vec,pch=19,col="red")
  points(1:42, m+s, col="blue", pch=19)
  
  abline(v=seq(7.5,42, by=7),col="green")
  
  plot(1:42,vec,pch=19,col="red", type = "l", main= paste("Day vs avg(mean) + period"), xlab= "day", ylab="counts")
  points(1:42, mean(m)+s, col="blue", pch=19, type="l")
  
  points(1:42,vec,pch=19,col="red")
  points(1:42, mean(m)+s, col="blue", pch=19)
  
  abline(v=seq(7.5,42, by=7),col="green")
}

plot_residuals <- function(resid, wd)
{
  x11()
  plot(1:42, resid, pch=19, col="red", type = "l", main= paste("Residuals"), xlab= "day", ylab="counts")
  points(1:42, resid, col="blue", pch=weather_data$rain+15)
  
  abline(v=seq(7.5,42, by=7),col="green")
  abline(h=0)
  
  legend(5, 4000, legend = c("Not rainy","Rainy"), pch=c(15,16))
}
##---------------------------------------------##

total_traffic <- data %>% group_by(NewDay) %>% summarize(n())
names(total_traffic) <- c("Day","N_travels")

## Compute mean and seasonality
m = mean_filter(total_traffic, q, d)
s = seasonal_filter(total_traffic, q, d, m)

## Plot mean and seasonality
plot_mean(total_traffic, m)
plot_mean_periodic(total_traffic, m, s)
##---------------------------------------------##

## Residuals
res = total_traffic$N_travels - m - s

## Plot residuals and rain
plot_residuals(res, weather_data)
# Rain tends to make residuals negative
# Sun teds to make residals positive

library(forecast)
x11()
autocorr = acf(res)
autocorr
## Apparently best WN or AR(1), no evident different model
## To few elements of the series
graphics.off()

#####---------------------------------------#####
#####           Hour-by-hour flow           #####
#####---------------------------------------#####
## Definition of datasets
sthour_data = data %>% group_by(startingt_hour, NewDay) %>% summarise(Flow = n()) %>% 
  ungroup() %>% rename(Str = startingt_hour, Day = NewDay)
head(sthour_data)

n_hours = length(unique(sthour_data$Str))
n_hours # 20

##---------------------------------------------##
## Plotting function
plot_hours <- function(stdata, weather_data, group)
{
  x11()
  
  poss_hours = unique(stdata$Str)
  
  if(group==FALSE)
  {
    par(mfrow=c(4,5))
    
    for(ii in c(1,6:24))
    {
      vec <- rep(0,42)
      for(jj in 1:42)
      {
        vec[jj] = as.integer(stdata %>% filter(Str==ii,Day==jj) %>% select(Flow))
        if(is.na(vec[jj]))
          vec[jj]=0
      }
      
      plot(1:42,vec,pch=19,col="red", type = "l", main= paste0("All stations together, h = ", ii), xlab= "day", ylab="counts")
      points(1:42, vec, col="blue", pch=weather_data$rain+15)
      abline(v=seq(7.5,42, by=7),col="green")
    }
  } else {
    par(mfrow = c(2,2))
    
    matx = matrix(rep(0,4*42),4,42)
    for(ii in c(1,6:24))
    {
      if(ii>=6&&ii<=9)
        kk=1
      if(ii>=10&&ii<=15)
        kk=2
      if(ii>=16&&ii<=19)
        kk=3
      if(ii>=20&&ii<=24 || ii==1)
        kk=4
      
      for(jj in 1:42)
      {
        val = as.integer(stdata %>% filter(Str==ii,Day==jj) %>% select(Flow))
        if(!is.na(val))
          matx[kk,jj] = matx[kk,jj]+val
      }
    }
    
    for (kk in 1:4)
    {
      plot(1:42,matx[kk,],pch=19,col="red", type = "l", main= paste0("All stations together, k = ", kk), xlab= "day", ylab="counts")
      points(1:42, matx[kk,], col="blue", pch=weather_data$rain+15)
      abline(v=seq(7.5,42, by=7),col="green")
    }
  }
}
##---------------------------------------------##

## Plotting per hour
plot_hours(sthour_data, weather_data, FALSE)

## Try to group the data in time slots
## {6,7,8,9}, {10,11,12,13,14,15}, {16,17,18,19}, {20,21,22,23,24,1}
plot_hours(sthour_data, weather_data, TRUE)
graphics.off()