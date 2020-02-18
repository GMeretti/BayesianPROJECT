############################################################
####### Poisson Base  ######################################
############################################################
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

Rt_class = read.csv('weat.csv')[,1]

########################################################################
################## Poisson all variables ##############################
########################################################################

Yyest <- embed(as.matrix(total_traffic)[,2],dimension = 8)[,2]
m1 <- c(0,1,1,1,1,0,1)
m1 <- diag(rep(m1,5))
Yyest_eq <- as.vector(m1%*%Yyest)
Yyest_diff <- Yyest - Yyest_eq
Y <- as.matrix(total_traffic)[8:42,2]
Yan <- embed(as.matrix(total_traffic)[,2],dimension = 8)[,8]
volume <- as.data.frame(cbind(Y,Yyest_eq,Yyest_diff,Yan, 
                              rep(c(rep(0,5),
                                    rep(1,2)),5),
                              embed(Rt_class,dimension=2)[7:41,], 
                              weather_data$Tmean[8:42]))
colnames(volume) <- c("Y","Yyest_eq","Yyest_diff","Yan","W","Rt","Ryest","T")

## Normalize old traffic for intepretatio of regression coefficients
volume[,2] <- volume[,2]/1000 
volume[,3] <- volume[,3] /1000
volume[,4] <- volume[,4] /1000

## Define the data 
dat = list(Y=volume$Y, X=volume[,2:8], n=35, r=7)

## A list of initial value for the MCMC algorithm in JAGS
#Classic
inits = function() {list(beta=rep(0,7), alpha = 0)}

## Define model
modelRegress=jags.model("Poisson.bug", data=dat, inits=inits,
                        n.adapt=10000, n.chains=3)

## Save model
save(modelRegress,file='Poisson_model.Rdata')

## Delete burnin
update(modelRegress, n.iter=19000) # this is the burn-in period

## tell JAGS to generate MCMC samples that we will use 
## to represent the posterior distribution 
variable.names <- c("beta", "alpha","Ypred", "lambda")

# their values recorded
n.iter=200000 
thin=20  

outputRegress = coda.samples(model=modelRegress, variable.names=variable.names,
                             n.iter=n.iter, thin=thin)
## Save the model
save(outputRegress,file= 'Poisson_output.Rdata')

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


alpha = as.data.frame(data.out[,grepl("alpha", names(data.out))])
x11()
par(mfrow=c(1,dim(alpha)[2]),mar=c(2.75,2.75,.5,.5),mgp=c(1.7,.7,0))
for (jj in 1:dim(alpha)[2])
{
  stationarity.plot(alpha[,jj],xlab="iteration",ylab=paste(expression(alpha),jj))
}

beta = as.data.frame(data.out[,grepl("beta", names(data.out))])
x11()
par(mfrow=c(1,dim(beta)[2]),mar=c(2.75,2.75,.5,.5),mgp=c(1.7,.7,0))
for (jj in 1:dim(beta)[2])
{
  stationarity.plot(beta[,jj],xlab="iteration",ylab=paste(expression(beta),jj))
}

lambda = as.data.frame(data.out[,grepl("lambda", names(data.out))])
lambda.plot = lambda[,c(1,5,10,15,20,25,30,35)]
x11()
par(mfrow=c(1,dim(lambda.plot)[2]),mar=c(2.75,2.75,.5,.5),mgp=c(1.7,.7,0))
for (jj in 1:dim(lambda.plot)[2])
{
  stationarity.plot(lambda.plot[,jj],xlab="iteration",ylab=paste(expression(lambda),jj))
}

#### 1 - check if beta parameters are significant, so if they are different from 0
beta_quantile <- apply(beta,2,quantile,c(0.05,0.95))
beta_m <- apply(beta,2,mean)
beta_quantile

# Plot for alpha
x11()
hist(alpha[,], nclass="fd", main="Intercept", xlab=expression(alpha), prob=T)
lines(density(alpha[,]),col="blue",lwd=2)
quantile(alpha[,],prob=c(0.025,0.5,0.975))
abline(v=quantile(alpha[,],prob=c(0.025)),col="red",lty=2,lwd=2)
abline(v=quantile(alpha[,],prob=c(0.5)),col="red",lty=1,lwd=2)
abline(v=quantile(alpha[,],prob=c(0.975)),col="red",lty=2,lwd=2)
legend("topright",legend=c("posterior median", "95% Credible bounds","kernel density smoother"),lwd=c(2,2,2), col=c("red","red","blue"),lty=c(1,2,1))

# Plot for beta
x11()
par(mfrow=c(3,3))
for (i in 1:7){
hist(beta[,i], nclass="fd", main=colnames(volume)[i+1], xlab=paste(expression(beta),i), prob=T)
lines(density(beta[,i]),col="blue",lwd=2)
quantile(beta[,i],prob=c(0.025,0.5,0.975))
abline(v=quantile(beta[,i],prob=c(0.025)),col="red",lty=2,lwd=2)
abline(v=quantile(beta[,i],prob=c(0.5)),col="red",lty=1,lwd=2)
abline(v=quantile(beta[,i],prob=c(0.975)),col="red",lty=2,lwd=2)
legend("topright",legend=c("posterior median", "95% Credible bounds","kernel density smoother"),lwd=c(2,2,2), col=c("red","red","blue"),lty=c(1,2,1))
}

################### predictions ##################################
pred_y <- as.data.frame(data.out[,grepl("Ypred", names(data.out))])
pred_int <- apply(pred_y,2,quantile,c(.05,.95))
pred_mean <- apply(pred_y,2,mean)
pred_sd <- apply(pred_y,2,sd)

## Plot of distribution of Y with posterior parameters, comparison with real Y - 90% credibility interval
sz = dim(data.out)[1]
x11()
plot(1:35, pred_mean, pch = 16, col = 'red', type ='o', xlab = 'day', ylab='traffic', main ='comparison between actual and predicted Y')
points(1:35, volume$Y, pch=15)
for(ii in 1:35)
{
  points(ii, quantile(pred_y[,ii], 0.05), col = 'darkred', pch = 16)
  points(ii, quantile(pred_y[,ii], 0.95), col = 'darkred', pch = 16)
  segments(ii, quantile(pred_y[,ii] , 0.05), ii, quantile(pred_y[,ii], 0.95), col  = 'darkred')
  if(ii!=35)
  {
    lines(c(ii, ii+1), c(quantile(pred_y[,ii] , 0.05), quantile(pred_y[,ii+1], 0.05)), col = 'darkred', lty = 3)
    lines(c(ii, ii+1), c(quantile(pred_y[,ii] , 0.95), quantile(pred_y[,ii+1], 0.95)), col = 'darkred', lty = 3)
  }
}



############# OUTLIERS ##########################
N_obs <- length(volume$Y)

outlier=(volume$Y>pred_int[2,]| volume$Y<pred_int[1,])
sum(outlier) # number of obs falling outside the credible intervals: 33 out of 35


########################################################################
################## Poisson Spike&Slab ##############################
########################################################################

## Define the data 
dat = list(Y=volume$Y, X=volume[,2:8], n=35, r=7)

## A list of initial value for the MCMC algorithm in JAGS
#SnS 1
inits = function() {list(gamma=rep(0,7), slab=rep(0,7), alpha = 0, theta=rep(0.5,7))}

## Define model
modelRegress=jags.model("Poisson_ss.bug", data=dat, inits=inits,
                        n.adapt=10000, n.chains=3)

## Save model
save(modelRegress,file='Poisson_ss_model.Rdata')

## Delete burnin
update(modelRegress, n.iter=19000) # this is the burn-in period

## tell JAGS to generate MCMC samples that we will use 
## to represent the posterior distribution 
variable.names <- c("beta", "alpha","Ypred", "lambda","gamma","theta")

# their values recorded
n.iter=200000 
thin=20  

outputRegress = coda.samples(model=modelRegress, variable.names=variable.names,
                             n.iter=n.iter, thin=thin)
## Save the model
save(outputRegress,file= 'Poisson_ss_output.Rdata')

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


alpha = as.data.frame(data.out[,grepl("alpha", names(data.out))])
x11()
par(mfrow=c(1,dim(alpha)[2]),mar=c(2.75,2.75,.5,.5),mgp=c(1.7,.7,0))
for (jj in 1:dim(alpha)[2])
{
  stationarity.plot(alpha[,jj],xlab="iteration",ylab=paste(expression(alpha),jj))
}

beta = as.data.frame(data.out[,grepl("beta", names(data.out))])
x11()
par(mfrow=c(1,dim(beta)[2]),mar=c(2.75,2.75,.5,.5),mgp=c(1.7,.7,0))
for (jj in 1:dim(beta)[2])
{
  stationarity.plot(beta[,jj],xlab="iteration",ylab=paste(expression(beta),jj))
}

gamma = as.data.frame(data.out[,grepl("gamma", names(data.out))])
x11()
par(mfrow=c(1,dim(gamma)[2]),mar=c(2.75,2.75,.5,.5),mgp=c(1.7,.7,0))
for (jj in 1:dim(gamma)[2])
{
  stationarity.plot(gamma[,jj],xlab="iteration",ylab=paste(expression(gamma),jj))
}

beta_quantile <- apply(beta,2,quantile,c(0.05,0.95))
beta_m <- apply(beta,2,mean)
beta_quantile


gamma_quantile <- apply(gamma,2,quantile,c(0.05,0.95))
gamma_m <- apply(gamma,2,mean)


# Plot for alpha
x11()
hist(alpha[,], nclass="fd", main="Intercept", xlab=expression(alpha), prob=T)
lines(density(alpha[,]),col="blue",lwd=2)
quantile(alpha[,],prob=c(0.025,0.5,0.975))
abline(v=quantile(alpha[,],prob=c(0.025)),col="red",lty=2,lwd=2)
abline(v=quantile(alpha[,],prob=c(0.5)),col="red",lty=1,lwd=2)
abline(v=quantile(alpha[,],prob=c(0.975)),col="red",lty=2,lwd=2)
legend("topright",legend=c("posterior median", "95% Credible bounds","kernel density smoother"),lwd=c(2,2,2), col=c("red","red","blue"),lty=c(1,2,1))

# Plot for beta
x11()
par(mfrow=c(3,3))
for (i in 1:7){
  hist(beta[,i], nclass="fd", main=colnames(volume)[i+1], xlab=paste(expression(beta),i), prob=T)
  lines(density(beta[,i]),col="blue",lwd=2)
  quantile(beta[,i],prob=c(0.025,0.5,0.975))
  abline(v=quantile(beta[,i],prob=c(0.025)),col="red",lty=2,lwd=2)
  abline(v=quantile(beta[,i],prob=c(0.5)),col="red",lty=1,lwd=2)
  abline(v=quantile(beta[,i],prob=c(0.975)),col="red",lty=2,lwd=2)
  legend("topright",legend=c("posterior median", "95% Credible bounds","kernel density smoother"),lwd=c(2,2,2), col=c("red","red","blue"),lty=c(1,2,1))
}

# MPM criterion
incl <- apply(gamma,2,mean)
incl

### -> remove Yan
#############################################
######### Without Yan #######################
#############################################
Yyest <- embed(as.matrix(total_traffic)[,2],dimension = 2)[,2]
m1 <- c(0,1,1,1,1,0,1)
m1 <- diag(c(c(1,1,1,1,0,1),rep(m1,5)))
Yyest_eq <- as.vector(m1%*%Yyest)
Yyest_diff <- Yyest - Yyest_eq
Y <- as.matrix(total_traffic)[2:42,2]
volume <- as.data.frame(cbind(Y,Yyest_eq,Yyest_diff, 
                              c(c(rep(0,4),rep(1,2)),rep(c(rep(0,5), rep(1,2)),5)),
                              embed(Rt_class,dimension=2)[1:41,], 
                              weather_data$Tmean[2:42]))
colnames(volume) <- c("Y","Yyest_eq","Yyest_diff","W","Rt","Ryest","T")

## Normalize old traffic for intepretatio of regression coefficients
volume[,2] <- volume[,2]/1000 
volume[,3] <- volume[,3] /1000

## Define the data 
dat = list(Y=volume$Y, X=volume[,2:7], n=41, r=6)
inits = function() {list(beta=rep(0,6), alpha = 0)}
## Define model
modelRegress=jags.model("Poisson.bug", data=dat, inits=inits,
                        n.adapt=10000, n.chains=3)
## Save model
save(modelRegress,file="Poisson_wYan_model.Rdata")

## Delete burnin
update(modelRegress, n.iter=19000) # this is the burn-in period

variable.names <- c("beta", "alpha","Ypred", "lambda")

# their values recorded
n.iter=200000 
thin=20  

outputRegress = coda.samples(model=modelRegress, variable.names=variable.names,
                             n.iter=n.iter, thin=thin)
## Save the model
save(outputRegress,file="Poisson_wYan_output.Rdata")

##### load Poisson_wYan_output -> Spike&Slab now is done, we have used it to remove Yan
load("Poisson_wYan_output.Rdata")

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


alpha = as.data.frame(data.out[,grepl("alpha", names(data.out))])
x11()
par(mfrow=c(1,dim(alpha)[2]),mar=c(2.75,2.75,.5,.5),mgp=c(1.7,.7,0))
for (jj in 1:dim(alpha)[2])
{
  stationarity.plot(alpha[,jj],xlab="iteration",ylab=paste(expression(alpha),jj))
}

beta = as.data.frame(data.out[,grepl("beta", names(data.out))])
x11()
par(mfrow=c(1,dim(beta)[2]),mar=c(2.75,2.75,.5,.5),mgp=c(1.7,.7,0))
for (jj in 1:dim(beta)[2])
{
  stationarity.plot(beta[,jj],xlab="iteration",ylab=paste(expression(beta),jj))
}

lambda = as.data.frame(data.out[,grepl("lambda", names(data.out))])
lambda.plot = lambda[,c(1,5,10,15,20,25,30,35)]
x11()
par(mfrow=c(1,dim(lambda.plot)[2]),mar=c(2.75,2.75,.5,.5),mgp=c(1.7,.7,0))
for (jj in 1:dim(lambda.plot)[2])
{
  stationarity.plot(lambda.plot[,jj],xlab="iteration",ylab=paste(expression(lambda),jj))
}


## Autocorrelation alpha
x11()
par(mfrow=c(1,dim(alpha)[2]))
for(i in 1:dim(alpha)[2]){
  acf(alpha[,i],lwd=3,col="red3", main=expression(alpha))
}

## Autocorrelation beta
x11()
par(mfrow=c(1,dim(beta)[2]))
for(i in 1:dim(beta)[2]){
  acf(beta[,i],lwd=3,col="red3", main=paste(expression(beta),i))
}


beta_quantile <- apply(beta,2,quantile,c(0.05,0.95))
beta_m <- apply(beta,2,mean)
beta_quantile

# Plot for alpha
x11()
hist(alpha[,], nclass="fd", main="Intercept", xlab=expression(alpha), prob=T)
lines(density(alpha[,]),col="blue",lwd=2)
quantile(alpha[,],prob=c(0.025,0.5,0.975))
abline(v=quantile(alpha[,],prob=c(0.025)),col="red",lty=2,lwd=2)
abline(v=quantile(alpha[,],prob=c(0.5)),col="red",lty=1,lwd=2)
abline(v=quantile(alpha[,],prob=c(0.975)),col="red",lty=2,lwd=2)
legend("topright",legend=c("posterior median", "95% Credible bounds","kernel density smoother"),lwd=c(2,2,2), col=c("red","red","blue"),lty=c(1,2,1))

# Plot for beta
nam <- colnames(volume)[c(2,3,5,6,7,8)]
x11()
par(mfrow=c(2,3))
for (i in 1:6){
  hist(beta[,i], nclass="fd", main=nam[i], xlab=paste(expression(beta),i), prob=T)
  lines(density(beta[,i]),col="blue",lwd=2)
  quantile(beta[,i],prob=c(0.025,0.5,0.975))
  abline(v=quantile(beta[,i],prob=c(0.025)),col="red",lty=2,lwd=2)
  abline(v=quantile(beta[,i],prob=c(0.5)),col="red",lty=1,lwd=2)
  abline(v=quantile(beta[,i],prob=c(0.975)),col="red",lty=2,lwd=2)
  legend("topright",legend=c("posterior median", "95% Credible bounds","kernel density smoother"),lwd=c(2,2,2), col=c("red","red","blue"),lty=c(1,2,1))
}


############## PREDICTION ##################################
pred_y <- as.data.frame(data.out[,grepl("Ypred", names(data.out))])
pred_int <- apply(pred_y,2,quantile,c(.05,.95))
pred_mean <- apply(pred_y,2,mean)
pred_sd <- apply(pred_y,2,sd)

## Plot of distribution of Y with posterior parameters, comparison with real Y - 90% credibility interval
sz = dim(data.out)[1]
x11()
plot(1:41, pred_mean, pch = 16, col = 'red', type ='o', xlab = 'day', ylab='traffic', main ='comparison between actual and predicted Y')
points(1:41, Y, pch=15)
for(ii in 1:41)
{
  points(ii, quantile(pred_y[,ii], 0.05), col = 'darkred', pch = 16)
  points(ii, quantile(pred_y[,ii], 0.95), col = 'darkred', pch = 16)
  segments(ii, quantile(pred_y[,ii] , 0.05), ii, quantile(pred_y[,ii], 0.95), col  = 'darkred')
  if(ii!=35)
  {
    lines(c(ii, ii+1), c(quantile(pred_y[,ii] , 0.05), quantile(pred_y[,ii+1], 0.05)), col = 'darkred', lty = 3)
    lines(c(ii, ii+1), c(quantile(pred_y[,ii] , 0.95), quantile(pred_y[,ii+1], 0.95)), col = 'darkred', lty = 3)
  }
}



######## OUTLIERS ############## 
N_obs <- length(volume$Y)

outlier=(volume$Y>pred_int[2,]| volume$Y<pred_int[1,])
sum(outlier) # number of obs falling outside the credible intervals: 33 out of 35

