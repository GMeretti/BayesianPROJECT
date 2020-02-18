###########################################################################
#####             Negative Binomial Model Mixed Effect   -  STAN      #####
###########################################################################
rm(list=ls())
library(coda)       
library(plotrix)
library(loo)
library(dplyr)
library(rjags)
library(rstan)
# for plots
library(ggplot2)
library(tidyr)
library(dplyr)
library(purrr)
library(ggsci)
require(gplots)
require(ggpubr)

# load dataset and weather, for computing the design matrix X
load("../../../Dataset/bikemi_data.RData")
data <- bike_nil_week_hour
head(data)

weather_data <- read.table("../../../Dataset/data_weather.csv", header = T, sep = ",")
head(weather_data)

total_traffic <- data %>% group_by(NewDay) %>% summarize(n())
names(total_traffic) <- c("Day","N_travels")

Y <- total_traffic$N_travels

X <- as.data.frame(cbind(rep(1,42), weather_data$rain, weather_data$Tmean))

G <- as.data.frame(cbind(rep(c(rep(1,5),0,0),6), rep(c(rep(0,5),1,1),6)))
colnames(G) <- c("WD","WE")

N <- length(Y)
p_fix <- 3
ngr <- 2

data_NB <-list(N = N, 
                p_fix = p_fix, 
                ngr = ngr,  
                Y = Y, 
                X = as.matrix(X), 
                G = as.matrix(G))

inits <- function() 
{
  list(beta = rep(0, p_fix), 
       theta = rep(0, ngr),
       sigma2_beta = rep(10, p_fix), 
       sigma2_theta = rep(10, ngr))
}

################### FIT THE MODEL

NEGBIN_GLMM <- stan(file = "NegBin_var.stan", 
                 data = data_NB,
                 chains = 2, 
                 iter = 30000, 
                 warmup = 10000, 
                 thin= 50, 
                 seed = 42, 
                 init = inits,
                 algorithm = 'NUTS')
save(NEGBIN_GLMM, file="NegBin_var_est.dat")
load("NegBin_var_est.dat")

rstan::traceplot(NEGBIN_GLMM, pars = "sigma2_beta", inc_warmup = TRUE)
rstan::traceplot(NEGBIN_GLMM, pars= "sigma2_theta", inc_warmup = TRUE)

rstan::traceplot(NEGBIN_GLMM, pars = "beta", inc_warmup = FALSE)
rstan::traceplot(NEGBIN_GLMM, pars= "theta", inc_warmup = FALSE)

coda_chain <- As.mcmc.list(NEGBIN_GLMM, pars = c("beta", "theta"))
summary(coda_chain)

# Posterior distributions -----------------------------------------------

# beta
plot_post <- NEGBIN_GLMM %>% 
  rstan::extract("beta") %>% 
  as.data.frame() %>% 
  map_df(as_data_frame, .id = 'param')

plot_post %>% 
  ggplot(aes(value, fill = param)) + 
  geom_density() + 
  facet_wrap(~param, scales = 'free') + 
  scale_fill_locuszoom() + 
  theme_minimal()+
  theme(legend.position="none")
 

# theta
plot_post <- NEGBIN_GLMM %>% 
  rstan::extract("theta") %>% 
  as.data.frame() %>% 
  map_df(as_data_frame, .id = 'param')

plot_post %>% 
  ggplot(aes(value, fill = param)) + 
  geom_density() + 
  facet_wrap(~param, scales = 'free') + 
  scale_fill_locuszoom() + 
  theme_minimal() + 
  coord_flip() +
  theme(legend.position="none")

############### PREDICTION ################################

params <- rstan::extract(NEGBIN_GLMM, c("beta", "theta", "r_param"), perm = T)
jpar <- cbind(params[[1]], params[[2]], params[[3]])
beta.data <- jpar[,c(1,2,3)]
theta.data <- jpar[,c(4,5)]
r_param.data <- jpar[,c(6,7)]
pred_y <- matrix(0,nrow=nrow(jpar),ncol = 42)

for(i in 1:42){
  for(j in 1:nrow(jpar)){
    mu <- exp(sum(beta.data[j,]*X[i,]) + sum(theta.data[j,]*G[i,]))
    r <- sum(r_param.data[j,]*G[i,])
    p <- r/(r+mu)
    pred_y[j,i] <- rnbinom(1,size = r, prob = p)
  }
}

pred_int <- apply(pred_y,2,quantile,c(.05,.95))
pred_mean <- apply(pred_y,2,mean)
pred_sd <- apply(pred_y,2,sd)

## Plot of distribution of Y with posterior parameters, comparison with real Y - 90% credibility interval
x11()
plot(1:42, pred_mean, pch = 16, ylim=c(0,20000), col = 'red', type ='o', xlab = 'day', ylab='traffic', main ='comparison between actual and predicted Y')
points(1:42, Y, pch=15)
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
N_obs <- length(Y)

outlier=(Y>pred_int[2,]| Y<pred_int[1,])
sum(outlier) # number of obs falling outside the credible intervals: 4 out of 35

