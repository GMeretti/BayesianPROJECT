###################################################################################
################# Weekend, Rain(1/0), Temperature: Negative binomial ##############
##################################################################################
rm(list=ls())
library(coda)       
library(plotrix)
library(loo)
library(dplyr)
library(rjags)
# load dataset and weather, for computing the design matrix X
load("../../../Dataset/bikemi_data.RData")
data <- bike_nil_week_hour
head(data)

weather_data <- read.table("../../../Dataset/data_weather.csv", header = T, sep = ",")
head(weather_data)

total_traffic <- data %>% group_by(NewDay) %>% summarize(n())
names(total_traffic) <- c("Day","N_travels")

Y <- as.matrix(total_traffic)[1:42,2]
volume <- as.data.frame(cbind(Y,rep(c(rep(0,5),rep(1,2)),6), weather_data$rain,weather_data$Tmean))
colnames(volume) <- c("Y","W","Rt","T")

## Define the data 
dat = list(Y=volume$Y, X=cbind(1,volume[,2:4]), n=42, mu.beta=rep(0,4), tau.beta=diag(rep(0.001,4)))

## A list of initial value for the MCMC algorithm in JAGS
inits = function() {list(beta=rep(0,4))}

## Define model
modelRegress=jags.model("NegBin.bug", data=dat, inits=inits,
                        n.adapt=100000, n.chains=3)

## Save model
save(modelRegress,file="NegBin_model.Rdata")

## Delete burnin
update(modelRegress, n.iter=200000) # this is the burn-in period

## tell JAGS to generate MCMC samples that we will use 
## to represent the posterior distribution 
variable.names <- c("beta", "Ypred", "p_param", "r_param", "lambda")

# their values recorded
n.iter=300000 
thin=100

outputRegress = coda.samples(model=modelRegress, variable.names=variable.names,
                             n.iter=n.iter, thin=thin)
## Save the model
save(outputRegress,file= 'NegBin_W_R_T_output.Rdata')

##---------------------------------------------##
#### Output analyisis
library(coda)       # per leggere catena mcmc
library(plotrix)    # per fare plot CIs

## Load trajectories
## the OUTPUT is mcmc.list object - coda-formatted object; 
## it needs to be converted into a matrix, in order to be "readable".
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


beta = as.data.frame(data.out[,grepl("beta", names(data.out))])
x11()
par(mfrow=c(1,dim(beta)[2]),mar=c(2.75,2.75,.5,.5),mgp=c(1.7,.7,0))
for (jj in 1:dim(beta)[2])
{
  stationarity.plot(beta[,jj],xlab="iteration",ylab=paste(expression(beta),jj))
}

r_param = as.data.frame(data.out[,grepl("r_param", names(data.out))])
x11()
par(mfrow=c(1,dim(r_param)[2]),mar=c(2.75,2.75,.5,.5),mgp=c(1.7,.7,0))
for (jj in 1:dim(r_param)[2])
{
  stationarity.plot(r_param[,jj],xlab="iteration",ylab=paste(expression(r),jj))
}

p_param = as.data.frame(data.out[,grepl("p_param", names(data.out))])
p_param.plot = p_param[,c(1,5,10,15,20,25,30,35)]
x11()
par(mfrow=c(1,dim(p_param.plot)[2]),mar=c(2.75,2.75,.5,.5),mgp=c(1.7,.7,0))
for (jj in 1:dim(p_param.plot)[2])
{
  stationarity.plot(p_param.plot[,jj],xlab="iteration",ylab=paste(expression(p),jj))
}

lambda = as.data.frame(data.out[,grepl("lambda", names(data.out))])
lambda.plot = lambda[,c(1,5,10,15,20,25,30,35)]
x11()
par(mfrow=c(1,dim(lambda.plot)[2]),mar=c(2.75,2.75,.5,.5),mgp=c(1.7,.7,0))
for (jj in 1:dim(lambda.plot)[2])
{
  stationarity.plot(lambda.plot[,jj],xlab="iteration",ylab=paste(expression(lambda),jj))
}

## Autocorrelation beta
x11()
par(mfrow=c(1,dim(beta)[2]))
for(i in 1:dim(beta)[2]){
  acf(beta[,i],lwd=3,col="red3", main=paste(expression(beta),i))
}

## Autocorrelation r
x11()
par(mfrow=c(1,dim(r_param)[2]))
for(i in 1:dim(r_param)[2]){
  acf(r_param[,i],lwd=3,col="red3", main=expression(r))
}


############### PARAMETERS ##################
beta_quantile <- apply(beta,2,quantile,c(0.05,0.95))
beta_m <- apply(beta,2,mean)
beta_quantile

# Plot for beta
nam <- colnames(volume)
nam[1] <- "Intercept"
x11()
par(mfrow=c(2,3))
for (i in 1:dim(beta)[2]){
  hist(beta[,i], nclass="fd", main=nam[i], xlab=paste(expression(beta),i), prob=T)
  lines(density(beta[,i]),col="blue",lwd=2)
  quantile(beta[,i],prob=c(0.025,0.5,0.975))
  abline(v=quantile(beta[,i],prob=c(0.025)),col="red",lty=2,lwd=2)
  abline(v=quantile(beta[,i],prob=c(0.5)),col="red",lty=1,lwd=2)
  abline(v=quantile(beta[,i],prob=c(0.975)),col="red",lty=2,lwd=2)
  legend("topright",legend=c("posterior median", "95% Credible bounds","kernel density smoother"),lwd=c(2,2,2), col=c("red","red","blue"),lty=c(1,2,1))
}


################## PREDICTIONS #########################À
pred_y <- as.data.frame(data.out[,grepl("Ypred", names(data.out))])
pred_int <- apply(pred_y,2,quantile,c(.05,.95))
pred_mean <- apply(pred_y,2,mean)
pred_sd <- apply(pred_y,2,sd)

## Plot of distribution of Y with posterior parameters, comparison with real Y - 90% credibility interval
sz = dim(data.out)[1]
x11()
plot(1:42, pred_mean, pch = 16, ylim=c(0,30000), col = 'red', type ='o', xlab = 'day', ylab='traffic', main ='comparison between actual and predicted Y')
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



####### OUTLIERS
N_obs <- length(volume$Y)

outlier=(volume$Y>pred_int[2,]| volume$Y<pred_int[1,])
sum(outlier) # number of obs falling outside the credible intervals: 2 out of 35
