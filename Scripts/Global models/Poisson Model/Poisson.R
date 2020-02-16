#################################################
#####             Poisson models            #####
#################################################
## Prerequisites
rm(list=ls())
library(rjags)   # to interface R with JAGS
library(pscl)
library(dplyr)

load("~/Desktop/Andrea/Bayes_prog/Bayesian_project/Dataset/bikemi_data.RData") 
data <- bike_nil_week_hour
rm(bike_nil_week_hour)
weather <- read.table("~/Desktop/Andrea/Bayes_prog/Bayesian_project/Dataset/data_weather.csv", header=T, sep=",") 

# extract total volume day by day
total_traffic <- data %>% group_by(NewDay) %>% summarize(n())
names(total_traffic) <- c("Day","N_travels")
##---------------------------------------------##

#################################################
#####        complete Poisson model         #####
#################################################
#Model name 
name <- "Poisson"

## Define companion dataset of the model
volume <- as.data.frame(cbind(embed(as.matrix(total_traffic)[,2],dimension = 8)[,c(1,2,8)], 
                              rep(c(rep(0,5),
                              rep(1,2)),5),
                              embed(weather$rain,dimension=2)[7:41,], 
                              weather$Tmean[8:42],
                              rep(c(rep(0,5),1,0),5),
                              rep(c(1,rep(0,6)),5)))
colnames(volume) <- c("Y","Yyest","Yan","W","Rt","Ryest","T","S","M")

## Normalize old traffic for intepretatio of regression coefficients
volume[,2] <- volume[,2]/1000 
volume[,3] <- volume[,3] /1000

## Define the data 
dat = list(Y=volume$Y, X=volume[,2:9], n=35, r=8)

## A list of initial value for the MCMC algorithm in JAGS
#Classic
inits = function() {list(beta=rep(0,8), alpha = 0)}

#SnS 1
inits = function() {list(gamma=rep(0,8), alpha = 0, delta=rep(0,8), prob=0.5)}

#SnS 2
inits = function() {list(beta=rep(0,8), alpha = 0, delta=rep(0,8), prob=0.5)}
## Define model
modelRegress=jags.model(paste0(name,".bug"), data=dat, inits=inits,
                        n.adapt=10000, n.chains=3)

## Save model
save(modelRegress,file=paste(name,'model.Rdata', sep = "_"))

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
save(outputRegress,file= paste(name,'output.Rdata', sep = "_"))

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

## Use either packege CODA, or standard tools in R
x11()
par(mfrow=c(3,3))
for(i in 1:9){
  acf(data.out[,paste0('beta.',i,'.')],lwd=3,
      col="red3",main=paste0("autocorrelation of beta",i))
}

#sub-chain containing the beta sample
beta.post <- data.out
#posterior mean of the beta parameters
beta.bayes  <- apply(beta.post,2,"mean")
beta.bayes

## Representation of the posterior chain of  beta0
chain <- beta.post[,3]
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
legend("topright",legend=c("posterior median", "95% Credible bounds","kernel density smoother"),
       lwd=c(2,2,2), col=c("red","red","blue"),lty=c(1,2,1))

graphics.off()
