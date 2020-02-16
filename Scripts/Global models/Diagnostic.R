#################################################
#####        BSTS models           #####
#################################################
## Prerequisites
rm(list=ls())
library(rjags)
library(dplyr)
library(bsts)
library(loo)
load("../../../Dataset/bikemi_data.RData")
data <- bike_nil_week_hour
head(data)

weather_data <- read.table("../../../Dataset/data_weather.csv", header = T, sep = ",")
head(weather_data)

total_traffic <- data %>% group_by(NewDay) %>% summarize(n())
names(total_traffic) <- c("Day","N_travels")

w = read.csv('weat.csv')
names(w) = c('index', 'mm', 'tmin', 'tmax')
##---------------------------------------------##

#################################################
#####         complete r-BSTS model         #####
#################################################
volume <- as.data.frame(cbind(total_traffic[,2], weather_data$rain, weather_data$Tmean))
colnames(volume) <- c("Y","Rt","T")
head(volume)



name <- "r_BSTS"
load(paste0(name,"_output.Rdata"))

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


alpha.data = as.data.frame(data.out[,grepl("alpha", names(data.out))])
x11()
par(mfrow=c(1,dim(alpha.data)[2]),mar=c(2.75,2.75,.5,.5),mgp=c(1.7,.7,0))
for (jj in 1:dim(alpha.data)[2])
{
  stationarity.plot(alpha.data[,jj],xlab="iteration",ylab=paste(expression(alpha),jj))
}

beta.data = as.data.frame(data.out[,grepl("beta", names(data.out))])
x11()
par(mfrow=c(1,dim(beta.data)[2]),mar=c(2.75,2.75,.5,.5),mgp=c(1.7,.7,0))
for (jj in 1:dim(beta.data)[2])
{
  stationarity.plot(beta.data[,jj],xlab="iteration",ylab=paste(expression(beta),jj))
}

nu.data = as.data.frame(data.out[,grepl("nu", names(data.out))])
x11()
par(mfrow=c(1,dim(nu.data)[2]),mar=c(2.75,2.75,.5,.5),mgp=c(1.7,.7,0))
for (jj in 1:dim(nu.data)[2])
{
  stationarity.plot(nu.data[,jj],xlab="iteration",ylab=paste(expression(nu),jj))
}

tau.data = as.data.frame(data.out[,grepl("tau", names(data.out))])
x11()
par(mfrow=c(1,dim(tau.data)[2]),mar=c(2.75,2.75,.5,.5),mgp=c(1.7,.7,0))
for (jj in 1:dim(tau.data)[2])
{
  stationarity.plot(tau.data[,jj],xlab="iteration",ylab=paste(expression(tau),jj))
}

rho.data = as.data.frame(data.out[,grepl("rho", names(data.out))])
rho.data = rho.data[,1:(dim(rho.data)[2]-1)]

mu.data = as.data.frame(data.out[,grepl("mu", names(data.out))])
mu.data = mu.data[,1:(dim(mu.data)[2]-1)]

delta.data = as.data.frame(data.out[,grepl("delta", names(data.out))])
delta.data = delta.data[,1:(dim(delta.data)[2]-1)]

gamma.data = as.data.frame(data.out[,grepl("gamma", names(data.out))])
gamma.data = gamma.data[,(1+lag):(dim(gamma.data)[2]-1)]

meanlik.data = as.data.frame(data.out[,grepl("meanlik", names(data.out))])

pred_y <- as.data.frame(data.out[,grepl("Ypred", names(data.out))]) # samples of the posterior for t = 1: 42
pred_int <- apply(pred_y,2,quantile,c(.05,.95))
pred_mean <- apply(pred_y,2,mean)
pred_sd <- apply(pred_y,2,sd)

graphics.off()

#######################
#######Error with predicetd expected values [mm_rain]
expected_values = colMeans(data.out)
expected.mu    = as.data.frame(expected_values[grepl("mu",    names(expected_values))])
expected.mu    = expected.mu[1:dim(expected.mu)[1]-1,]
expected.delta = as.data.frame(expected_values[grepl("delta", names(expected_values))])
expected.delta = expected.delta[1:dim(expected.delta)[1]-1,]
expected.gamma = as.data.frame(expected_values[grepl("gamma", names(expected_values))])
expected.gamma = expected.gamma[(lag+1):(dim(expected.gamma)[1]-1),]
expected.tau   = expected_values[grepl("tau",   names(expected_values))]
expected.beta  = expected_values[grepl("beta",  names(expected_values))]
expected.nu    = expected_values[grepl("nu",    names(expected_values))]
expected.alpha = expected_values[grepl("alpha", names(expected_values))]
pred = as.vector(expected.mu+expected.gamma+ (rep(c(rep(expected.beta[1],5),rep(expected.beta[3],2)),6))*f(w$mm) +(rep(c(rep(expected.beta[2],5),rep(expected.beta[4],2)),6))*volume$T) 
error = as.vector(pred- volume$Y)
error

## Plot of expected values[mm-rain]
x11()
plot(1:42, volume$Y, pch = 16, type ='o', xlab = 'day', ylab='traffic', main ='comparison between actual and predicted traffic')
points(1:42, expected.mu, col='blue', pch = 16)
lines(1:42, expected.mu, col = 'blue', lty=3)
points(1:42, expected.mu+expected.gamma, col='forestgreen', pch = 16)
lines(1:42, expected.mu+expected.gamma, col = 'forestgreen', lty=3)
points(1:42, pred, col='red', pch = 16)
lines(1:42, pred, col = 'red', lty=3)
########################

######################
## Error with predicetd expected values [no_temp]
expected_values = colMeans(data.out)
expected.mu    = as.data.frame(expected_values[grepl("mu",    names(expected_values))])
expected.mu    = expected.mu[1:dim(expected.mu)[1]-1,]
expected.delta = as.data.frame(expected_values[grepl("delta", names(expected_values))])
expected.delta = expected.delta[1:dim(expected.delta)[1]-1,]
expected.gamma = as.data.frame(expected_values[grepl("gamma", names(expected_values))])
expected.gamma = expected.gamma[(lag+1):(dim(expected.gamma)[1]-1),]
expected.rho   = as.data.frame(expected_values[grepl("rho",   names(expected_values))])
expected.rho   = expected.rho[1:dim(expected.rho)[1]-1,]
expected.tau   = expected_values[grepl("tau",   names(expected_values))]
expected.beta  = expected_values[grepl("beta",  names(expected_values))]
expected.nu    = expected_values[grepl("nu",    names(expected_values))]
expected.alpha = expected_values[grepl("alpha", names(expected_values))]
pred = as.vector(expected.mu+expected.gamma+expected.rho + (expected.beta[1]+expected.nu[1]*rep(c(rep(0,5),rep(1,2)),6))*w) 
error = as.vector(pred- volume$Y)
error

## Plot of expected values
x11()
plot(1:42, volume$Y, pch = 16, type ='o', xlab = 'day', ylab='traffic', main ='comparison between actual and predicted traffic')
points(1:42, expected.mu, col='blue', pch = 16)
lines(1:42, expected.mu, col = 'blue', lty=3)
points(1:42, expected.mu+expected.gamma, col='forestgreen', pch = 16)
lines(1:42, expected.mu+expected.gamma, col = 'forestgreen', lty=3)
points(1:42, expected.mu+expected.gamma+expected.rho, col='orange', pch = 16)
lines(1:42, expected.mu+expected.gamma+expected.rho, col = 'orange', lty=3)
points(1:42, pred, col='red', pch = 16)
lines(1:42, pred, col = 'red', lty=3)
## Plot of distribution of Y with posterior parameters, comparison with real Y - 90% credibility interval
sz = dim(data.out)[1]
x11()
plot(1:42, pred, pch = 16, col = 'red', type ='o', xlab = 'day', ylab='traffic', main ='comparison between actual and predicted Y')
points(1:42, volume$Y, pch=15)
for(ii in 1:42)
{
  points(ii, quantile(pred_y[,ii], 0.05), col = 'darkred', pch = 16)
  points(ii, quantile(pred_y[,ii], 0.95), col = 'darkred', pch = 16)
  segments(ii, quantile(pred_y[,ii] , 0.05), ii, quantile(pred_y[,ii], 0.95), col  = 'darkred')
  if(ii!=42)
  {
    lines(c(ii, ii+1), c(quantile(pred_y[,ii] , 0.05), quantile(pred_y[,ii+1], 0.05)), col = 'darkred', lty = 3)
    lines(c(ii, ii+1), c(quantile(pred_y[,ii] , 0.95), quantile(pred_y[,ii+1], 0.95)), col = 'darkred', lty = 3)
  }
}


## Representaton
## Representation of the posterior chain of  beta0
chain <- beta.data[,1]
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

N_obs <- length(volume$Y)
# A- OUTLIERS
outlier=(volume$Y>pred_int[2,]| volume$Y<pred_int[1,])
sum(outlier) # number of obs falling outside the credible intervals: 2 out of 42


# B - Bayesian predictive p-value
ind_pvalue=t((t(pred_y) > volume$Y))
pred_suptail=apply(ind_pvalue,2,sum)/n.chain
prob_tail<-rep(0, N_obs) 
for (i in 1 : N_obs){
  prob_tail[i]=min(pred_suptail[i],1-pred_suptail[i])
}

x11()
plot(1:42, prob_tail)
abline(h=0.05)
# OUTLIERS: those obs such that the BAYESIAN PREDICTIVE p-value 
#   is smaller than 0.05 (or than 0.1)
out2=(prob_tail<0.05)
sum(out2) # the same as before: 2 out of 42 outliers


# C - Bayesian residuals

bres <- (volume$Y - pred_mean)/pred_sd 
out_res = (abs(bres) > 2) #as a reference value we take 2

x11()
par(mfrow=c(1,1))
plot(1:42,bres,cex=1,ylim=c(-4,4))
abline(h=c(-2,2), lty = "dashed")

# E - WAIC
sigma <- 1/tau.data[,1]
loglik <- matrix(0,n.chain,N_obs)
for (i in 1:N_obs){
  loglik[,i] <- log(dnorm(volume$Y[i], meanlik.data[,i], sqrt(sigma)))
}
w <- waic(loglik)
w$estimates 

#or
loo(loglik)

# -INF results
# library(LaplacesDemon)
# LL <- pred_y
# WAIC(LL)
