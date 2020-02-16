#################################################
#####     SARIMA model for time series      #####
#################################################
## Prerequisites
rm(list=ls())
library(simts)
library(dplyr)
library(tseries)
load("../../Dataset/bikemi_data.RData")
data <- bike_nil_week_hour
set.seed(1996)
##---------------------------------------------##

time_series_volume_daybyday <- function()
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
  
  return(vec)
}

## Main plots
x11()
ts1 <- time_series_volume_daybyday()

volume_ts = gts(ts1, freq = 1, unit_ts = "travels", unit_time = "day", name_ts = "Volume of travels", data_name = "Evolution of travel volumes")
x11()
plot(volume_ts)
x11()
corr_analysis(volume_ts)
graphics.off()

## Differences
d_volume = gts(diff(volume_ts))
d2_volume = gts(diff(diff(volume_ts)))
d_s_volume = gts(diff(volume_ts, lag = 7))

x11()
par(mfrow = c(4,1))
plot(volume_ts)
plot(d_volume)
plot(d2_volume)
plot(d_s_volume)

x11()
corr_analysis(d_volume)
x11()
corr_analysis(d2_volume)
x11()
corr_analysis(d_s_volume)
graphics.off()
##---------------------------------------------##

## Model estimate
x11()
mod1 = estimate(SARIMA(ar = 2, i = 1, sar = 1, sma = 1, s = 7, si = 1), volume_ts, method = "mle")
check(mod1)

x11()
mod2 = estimate(SARIMA(ar = 2, i = 2, sar = 1, sma = 1, s = 7, si = 1), volume_ts, method = "mle")
check(mod2)

x11()
mod3 = estimate(SARIMA(ar = 1, i = 1, sar = 1, sma = 1, s = 7, si = 1), volume_ts, method = "mle")
check(mod3)

x11()
mod4 = estimate(SARIMA(ar=0, i=0,ma=-0.8, sar=0.95, si = 0 , sma = 0, s = 7), volume_ts, method = "mle") # sarma: arma(0,1)x(1,0)7 with predefined coeff
check(mod4)

x11()
mod5 = estimate(SARIMA(ar=0, i=0, ma=1, sar=1, si = 0 , sma = 0, s = 7), volume_ts, method = "mle") # sarma: arma(0,1)x(1,0)7
check(mod5)

x11()
mod6 = estimate(SARIMA(ar=0, i=0, ma=1, sar=1, si = 1, sma = 0, s = 7), volume_ts, method = "mle") # sarma: arma(0,1)x(1,1,0)7
check(mod6)

x11()
mod7 = estimate(SARIMA(ar=1, i=0, ma=0, sar=1, si = 0, sma = 0, s = 7), volume_ts, method = "mle") # sarma: arma(0,1)x(1,1,0)7
check(mod7)

graphics.off()

kpss.test(volume_ts, null="Trend")

predict(mod1, n.ahead = 5)
