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
#####         complete r-BSTS model         #####
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
var_gamma = diag(rep(0.00005, 6))

## Define the data 
dat = list(Y=volume$Y,  Rt=volume$Rt, T=volume$T, n=42, lag=5,
           mean_mu=mean_mu, mean_gamma=mean_gamma, var_gamma=var_gamma)

## A list of initial value for the MCMC algorithm in JAGS
inits = function() {list(beta=rep(0,2), alpha = 0)}

## Compute model
modelRegress=jags.model("r_BSTS.bug",data=dat,inits=inits,
                        n.adapt=10000,n.chains=2)

## Save the model
save(modelRegress, file='r_BSTS_model.Rdata') 

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
save(outputRegress,file='r_BSTS_output.Rdata')
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
par(mfrow=c(1,3),mar=c(2.75,2.75,.5,.5),mgp=c(1.7,.7,0))
stationarity.plot(data.out[,1],xlab="iteration",ylab=expression(alpha))
stationarity.plot(data.out[,2],xlab="iteration",ylab=expression(beta.1.))
stationarity.plot(data.out[,3],xlab="iteration",ylab=expression(beta.2.))


x11()
par(mfrow=c(1,5),mar=c(2.75,2.75,.5,.5),mgp=c(1.7,.7,0))
stationarity.plot(data.out[,181],xlab="iteration",ylab=expression(tau.1.))
stationarity.plot(data.out[,182],xlab="iteration",ylab=expression(tau.2.))
stationarity.plot(data.out[,183],xlab="iteration",ylab=expression(tau.3.))
stationarity.plot(data.out[,184],xlab="iteration",ylab=expression(tau.4.))
stationarity.plot(data.out[,185],xlab="iteration",ylab=expression(tau.5.))

graphics.off()

## Error with predicetd expected values
expected_values = colMeans(data.out)
pred = as.vector(expected_values[52:93] + expected_values[95:136]+ expected_values[138:179] + expected_values[2] * volume$Rt +expected_values[3]*volume$T) 
error = as.vector(pred- volume$Y)
error

## Plot of expected values
x11()
plot(1:42, volume$Y, pch = 16, type ='o', xlab = 'day', ylab='traffic', main ='comparison between actual and predicted traffic')
points(1:42, expected_values[95:136], col='blue', pch = 16)
lines(1:42, expected_values[95:136], col = 'blue', lty=3)
points(1:42, expected_values[95:136] + expected_values[52:93], col='forestgreen', pch = 16)
lines(1:42, expected_values[95:136] + expected_values[52:93], col = 'forestgreen', lty=3)
points(1:42, expected_values[95:136] + expected_values[52:93] + expected_values[138:179], col='orange', pch = 16)
lines(1:42, expected_values[95:136] + expected_values[52:93] + expected_values[138:179], col = 'orange', lty=3)
points(1:42, pred, col='red', pch = 16)
lines(1:42, pred, col = 'red', lty=3)

## Plot of uncertainty
x11()
par(mfrow  =c(2,2))
plot(1:42, expected_values[95:136], pch = 16, col = 'blue', type ='o', xlab = 'day', ylab='traffic', main ='comparison between actual and predicted mean')
for(ii in 1:42)
{
  points(ii, quantile(data.out[,94+ii], 0.05), col = 'dodgerblue', pch = 16)
  points(ii, quantile(data.out[,94+ii], 0.95), col = 'dodgerblue', pch = 16)
  segments(ii, quantile(data.out[,94+ii], 0.05), ii, quantile(data.out[,94+ii], 0.95), col  = 'dodgerblue')
  if(ii!=42)
  {
    lines(c(ii,ii+1), c(quantile(data.out[,94+ii], 0.05), quantile(data.out[,94+ii+1], 0.05)), col = 'dodgerblue', lty =3)
    lines(c(ii, ii+1), c(quantile(data.out[,94+ii], 0.95), quantile(data.out[,94+ii+1], 0.95)), col = 'dodgerblue', lty = 3)
  }
}

plot(1:42, expected_values[95:136] + expected_values[52:93], pch = 16, col = 'forestgreen', type ='o', xlab = 'day', ylab='traffic', main ='comparison between actual and predicted mean+period')
for(ii in 1:42)
{
  points(ii, quantile(data.out[,94+ii] + data.out[,51+ii], 0.05), col = 'greenyellow', pch = 16)
  points(ii, quantile(data.out[,94+ii] + data.out[,51+ii], 0.95), col = 'greenyellow', pch = 16)
  segments(ii, quantile(data.out[,94+ii] + data.out[,51+ii] , 0.05), ii, quantile(data.out[,94+ii]+ data.out[,51+ii], 0.95), col  = 'greenyellow')
  if(ii!=42)
  {
    lines(c(ii,ii+1), c(quantile(data.out[,94+ii] + data.out[,51+ii], 0.05), quantile(data.out[,94+ii+1] + data.out[,51+ii+1], 0.05)), col = 'greenyellow', lty =3)
    lines(c(ii, ii+1), c(quantile(data.out[,94+ii] + data.out[,51+ii], 0.95), quantile(data.out[,94+ii+1] + data.out[,51+ii+1], 0.95)), col = 'greenyellow', lty = 3)
  }
}

plot(1:42, expected_values[95:136] + expected_values[52:93] + expected_values[138:179], pch = 16, col = 'orange', type ='o', xlab = 'day', ylab='traffic', main ='comparison between actual and predicted mean+period+AR')
for(ii in 1:42)
{
  points(ii, quantile(data.out[,94+ii] + data.out[,51+ii] + data.out[,137+ii], 0.05), col = 'darkorange', pch = 16)
  points(ii, quantile(data.out[,94+ii] + data.out[,51+ii]+ data.out[,137+ii], 0.95), col = 'darkorange', pch = 16)
  segments(ii, quantile(data.out[,94+ii] + data.out[,51+ii]+ data.out[,137+ii] , 0.05), ii, quantile(data.out[,94+ii]+ data.out[,51+ii]+ data.out[,137+ii], 0.95), col  = 'darkorange')
  if(ii!=42)
  {
    lines(c(ii,ii+1), c(quantile(data.out[,94+ii] + data.out[,51+ii]+ data.out[,137+ii], 0.05), quantile(data.out[,94+ii+1] + data.out[,51+ii+1]+ data.out[,137+ii+1], 0.05)), col = 'darkorange', lty =3)
    lines(c(ii, ii+1), c(quantile(data.out[,94+ii] + data.out[,51+ii]+ data.out[,137+ii], 0.95), quantile(data.out[,94+ii+1] + data.out[,51+ii+1]+ data.out[,137+ii+1], 0.95)), col = 'darkorange', lty = 3)
  }
}

plot(1:42, pred, pch = 16, col = 'red', type ='o', xlab = 'day', ylab='traffic', main ='comparison between actual and predicted mean+period+AR+weather')
for(ii in 1:42)
{
  points(ii, quantile(data.out[,94+ii] + data.out[,51+ii] + data.out[,137+ii] + data.out[,2] * volume$Rt[ii] + data.out[,3]*volume$T[ii] , 0.05), col = 'darkred', pch = 16)
  points(ii, quantile(data.out[,94+ii] + data.out[,51+ii]+ data.out[,137+ii]+ data.out[,2] * volume$Rt[ii] + data.out[,3]*volume$T[ii] , 0.95), col = 'darkred', pch = 16)
  segments(ii, quantile(data.out[,94+ii] + data.out[,51+ii]+ data.out[,137+ii] + data.out[,2] * volume$Rt[ii] + data.out[,3]*volume$T[ii] , 0.05),
           ii, quantile(data.out[,94+ii]+ data.out[,51+ii]+ data.out[,137+ii] + data.out[,2] * volume$Rt[ii] + data.out[,3]*volume$T[ii] , 0.95), col  = 'darkred')
  if(ii!=42)
  {
    lines(c(ii,ii+1), c(quantile(data.out[,94+ii] + data.out[,51+ii]+ data.out[,137+ii] + data.out[,2] * volume$Rt[ii] + data.out[,3]*volume$T[ii] , 0.05),
                        quantile(data.out[,94+ii+1] + data.out[,51+ii+1]+ data.out[,137+ii+1] +  data.out[,2] * volume$Rt[ii+1] + data.out[,3]*volume$T[ii+1] , 0.05)), col = 'darkred', lty =3)
    lines(c(ii, ii+1), c(quantile(data.out[,94+ii] + data.out[,51+ii]+ data.out[,137+ii] + data.out[,2] * volume$Rt[ii] + data.out[,3]*volume$T[ii] , 0.95),
                         quantile(data.out[,94+ii+1] + data.out[,51+ii+1]+ data.out[,137+ii+1] + data.out[,2] * volume$Rt[ii+1] + data.out[,3]*volume$T[ii+1] , 0.95)), col = 'darkred', lty = 3)
  }
}

## Plot of distribution of Y with posterior parameters, comparison with real Y
sz = dim(data.out)[1]
with_noise = matrix(0, nrow = sz, ncol = 42)
n=10000

for(jj in 1:n)
{
  for(ii in 1:42)
  {
    with_noise[jj, ii] = data.out[jj,94+ii] + data.out[jj,51+ii] + data.out[jj,137+ii] + data.out[jj,2] * volume$Rt[ii] + data.out[jj,3]*volume$T[ii] + rnorm(n =1, mean = 0,  sd= 1/sqrt(data.out[jj,181]))
  }
}

future = c(43:49)
future_delta = matrix(0,n,8)
future_delta[,1] = data.out[1:n,46]
future_mu = matrix(0,n,8)
future_mu[,1] = data.out[1:n, 137]
future_gamma = matrix(0,n,13)
for(kk in 1:6)
{
future_gamma[, kk] = data.out[1:n,88+kk]
}
future_rho = matrix(0,n,8)
future_rho[,1] = data.out[1:n,180]
future_y = matrix(0,n,7)

datas = cbind(rbind(1,0,1,1,1,0,1), c(7,10,8,9,7,11,10))
for(ii in 1:n)
{
  for(jj in 1:7)
  {
    future_delta[ii, jj+1] = future_delta[ii, jj] + rnorm(n=1, mean=0, sd = 1/sqrt(data.out[ii, 183]))
    future_mu[ii, jj+1] = future_delta[ii,jj] + future_mu[ii,jj] + rnorm(n=1, mean=0, sd = 1/sqrt(data.out[ii, 182]))
    future_rho[ii,jj+1] = future_rho[ii,jj]*data.out[ii,1] + rnorm(n=1, mean= 0, sd = 1/sqrt(data.out[ii,185]))
    future_gamma[ii,jj+6] = -sum(future_gamma[ii,jj:(jj+5)]) + rnorm(n=1, mean=0, sd = 1/sqrt(data.out[ii,183]))
    future_y[ii,jj] = data.out[ii,2] * datas[jj,1] + data.out[ii,3]*datas[jj,2] + future_mu[ii,jj+1] + future_rho[ii,jj+1]+future_gamma[ii,jj+6]+ rnorm(n=1, mean=0, sd = 1/sqrt(data.out[ii,181]))
  }
}

x11()
plot(1:49, c(pred,colMeans(future_y)), pch = 16, col = 'red', type ='o', xlab = 'day', ylab='traffic', main ='comparison between actual and predicted Y')
points(1:42, volume$Y, pch=15)
for(ii in 1:42)
{
  points(ii, quantile(with_noise[1:n,ii], 0.05), col = 'darkred', pch = 16)
  points(ii, quantile(with_noise[1:n,ii], 0.95), col = 'darkred', pch = 16)
  segments(ii, quantile(with_noise[1:n,ii] , 0.05), ii, quantile(with_noise[1:n,ii], 0.95), col  = 'darkred')
  if(ii!=42)
  {
    lines(c(ii, ii+1), c(quantile(with_noise[1:n,ii] , 0.05), quantile(with_noise[1:n,ii+1], 0.05)), col = 'darkred', lty = 3)
    lines(c(ii, ii+1), c(quantile(with_noise[1:n,ii] , 0.95), quantile(with_noise[1:n,ii+1], 0.95)), col = 'darkred', lty = 3)
  }
}

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
