#Model with parameters node by node

rm(list=ls())
library(rjags)   # to interface R with JAGS
ins <- read.table("../../../Dataset/ins.csv", header = T, sep = ",")
outs <- read.table("../../../Dataset/outs.csv", header = T, sep = ",")
weather <- read.table("../../../Dataset/data_weather.csv", header=T, sep=",")


#The model 

p <- 3 #number of covariates + intercept
name <- "Poisson_nbn"

n <- dim(ins)[1]
m <- dim(ins)[2]
H <- n/42
week <- rep(c(rep(1,5), rep(2,2)), 6)

#the covariates
X <- matrix(c(rep(1,n), rep(weather$rain, each = H), rep(weather$Tmean, each = H), rep(1:H, 42), rep(week, each=H)), nrow = n, ncol = p+2)
#X <- matrix(c(rep(1,n), rep(1:H, 42), rep(week, each=H)), nrow = n, ncol = p+2)
I <- diag(rep(1,p))

dat = list(Nin=ins, Nout=outs, X=X, n=n, p=p, H=H, m = m, I = I, mu.0 = rep(0,p))
inits = function() {list(theta.0=rep(0,p), theta.1=rep(0,p), theta.2=rep(0,p),
                         alpha.0 = matrix(rep(0, 2*p), nrow = p, ncol=2),
                         alpha.1 = matrix(rep(0, 2*p), nrow = p, ncol=2),
                         alpha.2 = matrix(rep(0, 2*p), nrow = p, ncol=2),  
                         beta.1 = matrix(rep(0, H*p), nrow=p, ncol=H), beta.2 = matrix(rep(0, H*p), nrow=p, ncol=H),
                         beta.0 = matrix(rep(0, H*p), nrow=p, ncol=H),
                         phi = rep(0, m)
)}




#the modeling

modelRegress=jags.model(paste0(name,".bug"), data=dat, inits=inits,
                        n.adapt=10000, n.chains=3)

save(modelRegress,file=paste(name,'model.Rdata', sep = "_"))

## Delete burnin
update(modelRegress, n.iter=19000) # this is the burn-in period

## tell JAGS to generate MCMC samples that we will use 
## to represent the posterior distribution 
variable.names <- c("theta.0", "theta.1", "theta.2", "alpha.0","alpha.1","alpha.2","beta.0", "beta.1","beta.2", "phi", "X0_pred", "X1_pred", "X2_pred")

n.iter=10000 
thin=10  

outputRegress = coda.samples(model=modelRegress, variable.names=variable.names,
                             n.iter=n.iter, thin=thin)
## Save the model
save(outputRegress,file= paste(name,'output.Rdata', sep = "_"))


