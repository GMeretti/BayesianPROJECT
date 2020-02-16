#################################################
#####        residual BSTS models           #####
#################################################
## Prerequisites
rm(list=ls())
library(rjags)
library(dplyr)
library(bsts)
load("../../../Dataset/bikemi_data.RData")
data <- bike_nil_week_hour
head(data)

weather_data <- read.table("../../../Dataset/data_weather.csv", header = T, sep = ",")
head(weather_data)

total_traffic <- data %>% group_by(NewDay) %>% summarize(n())
names(total_traffic) <- c("Day","N_travels")
##---------------------------------------------##

#################################################
#####         T-less r-BSTS model          #####
#################################################
volume <- as.data.frame(cbind(total_traffic[,2], weather_data$rain, weather_data$Tmean))
colnames(volume) <- c("Y","Rt","T")
head(volume)

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

d = 7 # period is a week
q = 3

m = mean_filter(total_traffic, q, d)
s = seasonal_filter(total_traffic, q, d, m)

## Set initial conditons
mean_mu = mean(m)
mean_gamma = s[c(3:7,1)]
var_gamma = diag(rep(0.0001, 6))

## Define the data 
dat = list(Y=volume$Y,  Rt=volume$Rt, n=42, lag=5, S=7,
           mean_mu=mean_mu, mean_gamma=mean_gamma, var_gamma=var_gamma, dummy = 2*7)

## A list of initial value for the MCMC algorithm in JAGS
inits = function() {list(beta=0, alpha = 0)}

## Compute model
modelRegress=jags.model("Robust_T_less_r_BSTS.bug",data=dat,inits=inits,
                        n.adapt=10000,n.chains=2)

## Save the model
save(modelRegress, file='Robust_T_less_r_BSTS_model.Rdata') 

## Update burn-in
update(modelRegress,n.iter=20000) # this is the burn-in period

## MCMC samples to represent the posterior distribution 
variable.names <- c("beta", "alpha", "tau", "mu", "gamma", "rho", "delta")

# their values recorded
n.iter=500000
thin=20  

outputRegress=coda.samples(model=modelRegress,
                           variable.names=variable.names,n.iter=n.iter,thin=thin)

## Save output
save(outputRegress,file='Robust_T_less_r_BSTS_output.Rdata')
##-------------------------------------------##

#### Output analysis
library(coda)       
library(plotrix)

## Load trajectories as matrices
data.out=as.matrix(outputRegress) # trasform the mcmc.list into a matrix,
data.out=data.frame(data.out)     # or, better, into a dataframe (easiest to handle in R)
attach(data.out)
n.chain=dim(data.out)[1] 
n.chain # this is the final sample size

summary(data.out)
head(data.out)

## Validation of the stationarity
stationarity.plot<-function(x,...){
  S<-length(x)
  scan<-1:S
  ng<-min( round(S/100),10)
  group<-S*ceiling( ng*scan/S) /ng
  boxplot(x~group,...)               
}

x11()
par(mfrow=c(1,2),mar=c(2.75,2.75,.5,.5),mgp=c(1.7,.7,0))
stationarity.plot(data.out[,1],xlab="iteration",ylab=expression(alpha))
stationarity.plot(data.out[,2],xlab="iteration",ylab=expression(beta))


x11()
par(mfrow=c(1,5),mar=c(2.75,2.75,.5,.5),mgp=c(1.7,.7,0))
stationarity.plot(data.out[,208],xlab="iteration",ylab=expression(tau.1.))
stationarity.plot(data.out[,209],xlab="iteration",ylab=expression(tau.2.))
stationarity.plot(data.out[,210],xlab="iteration",ylab=expression(tau.3.))
stationarity.plot(data.out[,211],xlab="iteration",ylab=expression(tau.4.))
stationarity.plot(data.out[,212],xlab="iteration",ylab=expression(tau.5.))

graphics.off()

## Error with predicetd expected values
expected_values = colMeans(data.out)
pred = as.vector(expected_values[122:163] + expected_values[65:106]+ expected_values[165:206] + expected_values[2] * volume$Rt) 
error = as.vector(pred- volume$Y)
error

## Plot of expected values
x11()
plot(1:42, volume$Y, pch = 16, type ='o', xlab = 'day', ylab='traffic', main ='comparison between actual and predicted traffic')
points(1:42, expected_values[122:163], col='blue', pch = 16)
lines(1:42, expected_values[122:163], col = 'blue', lty=3)
points(1:42, expected_values[122:163] + expected_values[65:106], col='forestgreen', pch = 16)
lines(1:42, expected_values[122:163] + expected_values[65:106], col = 'forestgreen', lty=3)
points(1:42, expected_values[122:163] + expected_values[65:106] + expected_values[165:206], col='orange', pch = 16)
lines(1:42, expected_values[122:163] + expected_values[65:106] + expected_values[165:206], col = 'orange', lty=3)
points(1:42, pred, col='red', pch = 16)
lines(1:42, pred, col = 'red', lty=3)


## Plot of distribution of Y with posterior parameters, comparison with real Y
sz = dim(data.out)[1]
with_noise = matrix(0, nrow = sz, ncol = 42)
n=10000
S=7

for(jj in 1:n)
{
  for(ii in 1:42)
  {
    with_noise[jj, ii] = data.out[jj,121+ii] + data.out[jj,64+ii] + data.out[jj,164+ii] + data.out[jj,2] * volume$Rt[ii] + rnorm(n =1, mean = 0,  sd= 1/sqrt(data.out[jj,208]))
  }
}

future_delta = matrix(0,n,3*S)
future_mu = matrix(0,n,3*S)
for (kk in 1:(2*S))
{
  future_delta[,kk] = data.out[1:n, 44+kk]
  future_mu[,kk] = data.out[1:n, 149+kk]
}
future_gamma = matrix(0,n,2*S-1)
for(kk in 1:6)
{
  future_gamma[, kk] = data.out[1:n, 100+kk]
}

future_rho = matrix(0,n,S+1)
future_rho[,1] = data.out[1:n, 206]

future = 43:49
future_y = matrix(0,n,7)

datas = cbind(rbind(1,0,1,1,1,0,1))
for(ii in 1:n)
{
  for(jj in 1:S)
  {
    future_delta[ii, jj+2*S] = 1/2*future_delta[ii, jj+2*S-1] + 1/3*future_delta[ii, jj+S] + 1/6*future_delta[ii, jj] + rnorm(n=1, mean=0, sd = 1/sqrt(data.out[ii, 210]))
    future_mu[ii, jj+2*S] = 1/2*future_delta[ii,jj+2*S-1] + 1/2*future_mu[ii,jj+2*S-1] +1/3*future_delta[ii,jj+S] + 1/3*future_mu[ii,jj+S]+1/6*future_delta[ii,jj] + 1/6*future_mu[ii,jj] + rnorm(n=1, mean=0, sd = 1/sqrt(data.out[ii, 209]))
    future_rho[ii,jj+1] = future_rho[ii,jj]*data.out[ii,1] + rnorm(n=1, mean= 0, sd = 1/sqrt(data.out[ii,212]))
    future_gamma[ii,jj+S-1] = -sum(future_gamma[ii,jj:(jj+S-2)]) + rnorm(n=1, mean=0, sd = 1/sqrt(data.out[ii,211]))
    future_y[ii,jj] = data.out[ii,2] * datas[jj,1] + future_mu[ii,jj+2*S] + future_rho[ii,jj+1]+future_gamma[ii,jj+S-1]+ rnorm(n=1, mean=0, sd = 1/sqrt(data.out[ii,208]))
  }
}

x11()
plot(1:49, c(pred,colMeans(future_y)), pch = 16, col = 'red', type ='o', xlab = 'day', ylab='traffic', main ='comparison between actual and predicted Y')
points(1:42, volume$Y, pch=15)

with_noise_2 = cbind(with_noise[1:n,], future_y[1:n,])
for(ii in 1:49)
{
  points(ii, quantile(with_noise_2[1:n,ii], 0.05), col = 'darkred', pch = 16)
  points(ii, quantile(with_noise_2[1:n,ii], 0.95), col = 'darkred', pch = 16)
  segments(ii, quantile(with_noise_2[1:n,ii] , 0.05), ii, quantile(with_noise_2[1:n,ii], 0.95), col  = 'darkred')
  if(ii!=49)
  {
    lines(c(ii, ii+1), c(quantile(with_noise_2[1:n,ii] , 0.05), quantile(with_noise_2[1:n,ii+1], 0.05)), col = 'darkred', lty = 3)
    lines(c(ii, ii+1), c(quantile(with_noise_2[1:n,ii] , 0.95), quantile(with_noise_2[1:n,ii+1], 0.95)), col = 'darkred', lty = 3)
  }
}

## LORENZO TODO
## Plots
x11()
par(mfrow=c(3,3))
for(i in 1:8){
  acf(data.out[,i],lwd=3,col="red3")
}

## Representaton
## Representation of the posterior chain of  beta0
chain <- data.out[,2]
#Divide the plot device in three sub-graph regions
#two square on the upper and a rectangle on the bottom
x11()
layout(matrix(c(1,2,3,3),2,2,byrow=T))
#trace-plot of the posterior chain
plot(chain,type="l",main="Trace plot")
# autocorrelation plot
acf(chain,lwd=3,col="red3",main="autocorrelation")
#Histogram
hist(chain,nclass="fd",freq=F,main="Posterior",col="gray") 
## Overlap the kernel-density 
lines(density(chain),col="blue",lwd=2)
## Posterior credible interval of beta0
quantile(chain,prob=c(0.025,0.5,0.975))

## Display the posterior credible interval on the graph
abline(v=quantile(chain,prob=c(0.025)),col="red",lty=2,lwd=2)
abline(v=quantile(chain,prob=c(0.5)),col="red",lty=1,lwd=2)
abline(v=quantile(chain,prob=c(0.975)),col="red",lty=2,lwd=2)
## Add the legend to the plot
legend("topright",legend=c("posterior median", "95% Credible bounds","kernel density smoother"),lwd=c(2,2,2), col=c("red","red","blue"),lty=c(1,2,1))

graphics.off()
