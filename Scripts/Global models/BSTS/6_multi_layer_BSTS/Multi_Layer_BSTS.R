#################################################
#####        residual BSTS models           #####
#################################################
## Prerequisites
rm(list=ls())
library(rjags)
library(dplyr)
library(bsts)
library(loo)
load("../../../../Dataset/bikemi_data.RData")
data <- bike_nil_week_hour
head(data)

weather_data <- read.table("../../../../Dataset/data_weather.csv", header = T, sep = ",")
head(weather_data)

total_traffic <- data %>% group_by(NewDay) %>% summarize(n())
names(total_traffic) <- c("Day","N_travels")

w = read.csv('weat.csv')[,1]
##---------------------------------------------##

#################################################
#####          Robust r-BSTS model          #####
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

get_hours <- function(stdata)
{
  poss_hours = unique(stdata$Str)
  par(mfrow = c(2,2))
  
  matx = matrix(rep(0,4*42),4,42)
  for(ii in c(1, 6:24))
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
      is_any = dim(stdata %>% filter(Str==ii,Day==jj))[1]
      if(is_any!=0)
      {
        val = as.integer((stdata %>% filter(Str==ii,Day==jj))[3])
        matx[kk,jj] = matx[kk,jj]+val
      }
    }
  }
  
  return(matx)
}

d = 7 # period is a week
q = 3

m = mean_filter(total_traffic, q, d)
s = seasonal_filter(total_traffic, q, d, m)

sthour_data = data %>% group_by(startingt_hour, NewDay) %>% summarise(Flow = n()) %>% 
  ungroup() %>% rename(Str = startingt_hour, Day = NewDay)
head(sthour_data)

n_hours = length(unique(sthour_data$Str))
n_hours # 20

matr = get_hours(sthour_data)

Y = c()
for(jj in 1:42)
{
  Y = c(Y,matr[,jj])
}
Y

v = rep(0,4)
count = 0
for(jj in 1:42)
{
  for(kk in 1:4)
  {
    if(jj%%7!=0 && jj%%7!=6)
    {
      count = count+1
      v[kk] = v[kk] + matr[kk, jj] - (m[jj]+s[jj])/4
    }
  }
}
v = v/(count/4)
  
## Set initial conditons
mean_mu = mean(m)/4
mean_gamma = s[c(3:7,1)]/4
var_gamma = diag(rep(0.0001, 6))
var_chi = diag(rep(0.001,3))
mean_chi = c(v[3:4],v[1])

## Define the data 
dat = list(Y=Y,  Rt=w, T=volume$T, dummy_weeks=2, tot_weeks=6, S=7, F=4,
           mean_mu=mean_mu, mean_gamma=mean_gamma, Wd = c(rep(0,5),rep(1,2)),
           mean_chi=mean_chi)

## A list of initial value for the MCMC algorithm in JAGS
inits = function() {list(beta=rep(0,2))}

## Compute model
modelRegress=jags.model("Multi_Layer_BSTS.bug",data=dat,inits=inits,
                        n.adapt=10000,n.chains=2)

## Save the model
save(modelRegress, file='Robust_Multi_Layer_r_BSTS_model.Rdata') 

## Update burn-in
update(modelRegress,n.iter=20000) # this is the burn-in period

## MCMC samples to represent the posterior distribution 
variable.names <- c("beta", "tau", "mu", "gamma", "delta", "chi", 'xi', "nu","Ypred","meanlik")

# their values recorded
n.iter=1000000
thin=100 

outputRegress=coda.samples(model=modelRegress,
                           variable.names=variable.names,n.iter=n.iter,thin=thin)

## Save output
save(outputRegress,file='Robust_Multi_Layer_r_BSTS_output.Rdata')
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

xi.data = as.data.frame(data.out[,grepl("xi", names(data.out))])
x11()
par(mfrow=c(1,dim(xi.data)[2]),mar=c(2.75,2.75,.5,.5),mgp=c(1.7,.7,0))
for (jj in 1:dim(xi.data)[2])
{
  stationarity.plot(xi.data[,jj],xlab="iteration",ylab=paste(expression(xi),jj))
}


tau.data = as.data.frame(data.out[,grepl("tau", names(data.out))])
x11()
par(mfrow=c(1,dim(tau.data)[2]),mar=c(2.75,2.75,.5,.5),mgp=c(1.7,.7,0))
for (jj in 1:dim(tau.data)[2])
{
  stationarity.plot(tau.data[,jj],xlab="iteration",ylab=paste(expression(tau),jj))
}

graphics.off()

## Autocorrelation beta
x11()
par(mfrow=c(1,dim(beta.data)[2]))
for(i in 1:dim(beta.data)[2]){
  acf(beta.data[,i],lwd=3,col="red3", main=paste(expression(beta),i))
}

## Autocorrelation nu
x11()
par(mfrow=c(1,dim(nu.data)[2]))
for(i in 1:dim(nu.data)[2]){
  acf(nu.data[,i],lwd=3,col="red3", main=paste(expression(nu),i))
}

## Autocorrelation tau
x11()
par(mfrow=c(1,dim(tau.data)[2]))
for(i in 1:dim(tau.data)[2]){
  acf(tau.data[,i],lwd=3,col="red3", main=paste(expression(tau),i))
}

graphics.off()

mu.data = as.data.frame(data.out[,grepl("mu", names(data.out))])
mu.data = mu.data[,(dim(mu.data)[2]-42):(dim(mu.data)[2]-1)]
delta.data = as.data.frame(data.out[,grepl("delta", names(data.out))])
delta.data = delta.data[,(dim(delta.data)[2]-42):(dim(delta.data)[2]-1)]
gamma.data = as.data.frame(data.out[,grepl("gamma", names(data.out))])
gamma.data = gamma.data[,(dim(gamma.data)[2]-42):(dim(gamma.data)[2]-1)]
chi.data = as.data.frame(data.out[,grepl("chi", names(data.out))])
chi.data = chi.data[,(dim(chi.data)[2]-168):(dim(chi.data)[2]-1)]

meanlik.data = as.data.frame(data.out[,grepl("meanlik", names(data.out))])

pred_y <- as.data.frame(data.out[,grepl("Ypred", names(data.out))]) # samples of the posterior for t = 1: 42
pred_int <- apply(pred_y,2,quantile,c(.05,.95))
pred_mean <- apply(pred_y,2,mean)
pred_sd <- apply(pred_y,2,sd)

## Error with predicetd expected values DEPRECATED
expected_values = colMeans(data.out)
expected.mu    = as.data.frame(expected_values[grepl("mu",    names(expected_values))])
expected.mu    = expected.mu[(dim(expected.mu)[1]-42):(dim(expected.mu)[1]-1),]
expected.mu    = rep(expected.mu, each=4)
expected.delta = as.data.frame(expected_values[grepl("delta", names(expected_values))])
expected.delta = expected.delta[(dim(expected.delta)[1]-42):(dim(expected.delta)[1]-1),]
expected.delta = rep(expected.delta, each=4)
expected.gamma = as.data.frame(expected_values[grepl("gamma", names(expected_values))])
expected.gamma = expected.gamma[(dim(expected.gamma)[1]-42):(dim(expected.gamma)[1]-1),]
expected.gamma = rep(expected.gamma, each=4)
expected.chi   = as.data.frame(expected_values[grepl("chi", names(expected_values))])
expected.chi   = expected.chi[(dim(expected.chi)[1]-168):(dim(expected.chi)[1]-1),]
expected.xi    = as.data.frame(expected_values[grepl("xi", names(expected_values))])
expected.xi    = expected.xi[1:4,]
expected.xi    = rep(expected.xi,42)
expected.tau   = expected_values[grepl("tau",   names(expected_values))]
expected.beta  = expected_values[grepl("beta",  names(expected_values))]
expected.nu    = expected_values[grepl("nu",    names(expected_values))]
expected.coef1 = rep(expected.beta[1],168)+rep(rep(c(rep(0,5),rep(expected.nu[1],2)),6),each=4) + expected.xi
expected.coef2 = rep(expected.beta[2],168)+rep(rep(c(rep(0,5),rep(expected.nu[2],2)),6),each=4)
pred =  as.vector(expected.mu+ expected.gamma + expected.chi + expected.coef1*rep(w,each=4) + expected.coef2*rep(volume$T,each=4))
error = as.vector(pred- Y)
error

## Plot of distribution of Y with posterior parameters, comparison with real Y - 90% credibility interval
## Use this nstead
n=1000
sz = dim(data.out)[1]
with_noise = matrix(0, nrow = n, ncol = 168)


for(jj in 1:n)
{
  for(ii in 1:42)
  {
    for(kk in 1:4)
    {
      with_noise[jj, (ii-1)*4+kk] = max(mu.data[jj,ii] + gamma.data[jj,ii] + chi.data[jj,kk+(ii-1)*4]+ (beta.data[jj,1]+nu.data[jj,1]*((ii%%7)==0||(ii%%7)==6)+xi.data[jj,kk])*w[ii] + (beta.data[jj,2]+nu.data[jj,2]*((ii%%7)==0||(ii%%7)==6))*volume$T[ii] + rnorm(n =1, mean = 0,  sd= 1/sqrt(tau.data[jj,1])),0)
    }
  }
}

x11()
plot(1:168, colMeans(with_noise), pch = 16, col = 'red', type ='o', xlab = 'day', ylab='traffic', main ='comparison between actual and predicted Y', lwd=2)
points(1:168, Y, pch=15)
lines(1:168, Y, lwd=2)

for(ii in 1:168)
{
  points(ii, quantile(with_noise[1:n,ii], 0.05), col = 'darkred', pch = 16)
  points(ii, quantile(with_noise[1:n,ii], 0.95), col = 'darkred', pch = 16)
  segments(ii, quantile(with_noise[1:n,ii] , 0.05), ii, quantile(with_noise[1:n,ii], 0.95), col  = 'darkred', lty=2)
  if(ii%%4==0)
  {
    abline(v=ii+0.5, pch=0.5, col='brown')
  }
}

#Evaluation
loo(loglik)


## Representaton
## Representation of the posterior chain of  beta0
chain <- mu.data[,2]
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