############################################################
####### Predictive GOF analysis: DAY-BY-DAY TRAVELS ########
############################################################
library(coda)       
library(plotrix)
# load dataset and weather, for computing the design matrix X
load("../../Dataset/bikemi_data.RData")
data <- bike_nil_week_hour
head(data)

weather_data <- read.table("../../Dataset/data_weather.csv", header = T, sep = ",")
head(weather_data)

total_traffic <- data %>% group_by(NewDay) %>% summarize(n())
names(total_traffic) <- c("Day","N_travels")

volume <- as.data.frame(cbind(total_traffic[,2], weather_data$rain, weather_data$Tmean))
colnames(volume) <- c("Y","Rt","T")
head(volume)

#################################################
#####         complete r-BSTS model         #####
#################################################

load("r-BSTS Model/r_BSTS_output.Rdata")
## Load trajectories as matrices
data.out=as.matrix(outputRegress) # trasform the mcmc.list into a matrix,
data.out=data.frame(data.out)     # or, better, into a dataframe (easiest to handle in R)
attach(data.out)
n.chain=dim(data.out)[1] 
n.chain # this is the final sample size

summary(data.out)
head(data.out)


#### 1 - check if beta parameters are significant, so if they are different from 0
# Plot for beta1
x11()
hist(beta.1., nclass="fd", main="RAIN", xlab=expression(beta[1]), prob=T)
lines(density(beta.1.),col="blue",lwd=2)
quantile(beta.1.,prob=c(0.025,0.5,0.975))
abline(v=quantile(beta.1.,prob=c(0.025)),col="red",lty=2,lwd=2)
abline(v=quantile(beta.1.,prob=c(0.5)),col="red",lty=1,lwd=2)
abline(v=quantile(beta.1.,prob=c(0.975)),col="red",lty=2,lwd=2)
legend("topright",legend=c("posterior median", "95% Credible bounds","kernel density smoother"),lwd=c(2,2,2), col=c("red","red","blue"),lty=c(1,2,1))
# Plot for beta2
x11()
hist(beta.2., nclass="fd", main="TEMPERATURE", xlab=expression(beta[2]), prob=T)
lines(density(beta.2.),col="blue",lwd=2)
quantile(beta.2.,prob=c(0.025,0.5,0.975))
abline(v=quantile(beta.2.,prob=c(0.025)),col="red",lty=2,lwd=2)
abline(v=quantile(beta.2.,prob=c(0.5)),col="red",lty=1,lwd=2)
abline(v=quantile(beta.2.,prob=c(0.975)),col="red",lty=2,lwd=2)
legend("topright",legend=c("posterior median", "95% Credible bounds","kernel density smoother"),lwd=c(2,2,2), col=c("red","red","blue"),lty=c(1,2,1))

##### 2 - evolution of the trend 
trend <- data.out[,138:179]
m_trend <- apply(trend,2,mean)
q_trend <- apply(trend,2,quantile,c(0.05,0.95))
# Plot the 90% CI and the datapoints
x11()
par(mfrow=c(1,1))
ind <- 1:length(volume$Y)
matplot(rbind(ind,ind), q_trend,type="l",lty=1,col=1,xlab="day",ylab=expression(mu))
lines(ind, m_trend, col="red")
abline(v=seq(7.5,42.5, by=7), col="blue")

#### 3 - some plot of the posterior distributions of y(t)s (for t = 1, 2 for example)
volume$Y # given obesrvations for t = 1: 42 
pred_y <- data.out[,1:42] # samples of the posterior for t = 1: 42
pred_int <- apply(pred_y,2,quantile,c(.05,.95))
pred_mean <- apply(pred_y,2,mean)
pred_sd <- apply(pred_y,2,sd)

x11()
par(mfrow=c(2,1))
hist(pred_y[,1], nclass="fd", main="Y1", prob=T)
points(volume$Y[1],0, col="red")
abline(v=pred_int[,1], col="green")
abline(v=pred_mean[1], lty="dashed", col="blue")

hist(pred_y[,2], nclass="fd", main="Y2", prob=T)
points(volume$Y[2],0, col="red")
abline(v=pred_int[,2], col="green")
abline(v=pred_mean[2], lty="dashed", col="blue")


####### 3 - PREDICTIVE GOODNESS-OF-FIT 
N_obs <- length(volume$Y)

# A - Plot the 90% CI and the datapoints
x11()
par(mfrow=c(1,1))
ind <- 1:N_obs
matplot(rbind(ind,ind), pred_int, type="l", lty=1, col=1, xlab="day", ylab="travels")
points(ind, volume$Y, pch=19, cex=0.2, col="red")
# OUTLIERS
outlier=(volume$Y>pred_int[2,]| volume$Y<pred_int[1,])
sum(outlier) # number of obs falling outside the credible intervals: 36 out of 42

# B - Bayesian predictive p-value
ind_pvalue=t((t(pred_y) > volume$Y))
pred_suptail=apply(ind_pvalue,2,sum)/n.chain
prob_tail<-rep(0, N_obs) 
for (i in 1 : N_obs){
  prob_tail[i]=min(pred_suptail[i],1-pred_suptail[i])
}

x11()
plot(ind, prob_tail)
abline(h=0.05)

# OUTLIERS: those obs such that the BAYESIAN PREDICTIVE p-value 
#   is smaller than 0.05 (or than 0.1)
out2=(prob_tail<0.05)
sum(out2) # the same as before: 36 out of 42 outliers

# C - Bayesian residuals

bres <- (volume$Y - pred_mean)/pred_sd 
out_res = (abs(bres) > 2) #as a reference value we take 2

x11()
par(mfrow=c(1,1))
plot(ind,bres,cex=1)
abline(h=c(-2,2), lty = "dashed")


# D - CPO (sort of cross-validation)
alp <- data.out$alpha
delta <- data.out[,46:88]
gamma <- data.out[,89:136]
mu <- data.out[,137:179]
rho <- data.out[,180:222]
tau <- data.out[,223:227]
sigma <- 1 / tau[,1] # tau[,1] is tau(eps)

cpo <- seq(1,N_obs)
for (i in 1:N_obs){
  m <- rho[,i+1] + mu[,i+1] + gamma[1,i+6] + beta.1.*volume$Rt[i] + beta.2.*volume$T[i]
  cpo[i]=1/mean(1/dnorm(volume$Y[i], m, sqrt(sigma)))
}

LPML <- sum(log(cpo))
LPML

