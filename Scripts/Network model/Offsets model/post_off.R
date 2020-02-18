#Analysis of the posterior

name <- "Poisson_offsets"

load(paste(name, "output.Rdata", sep = "_"))

ins <- read.table("../../../Dataset/ins.csv", header = T, sep = ",")
outs <- read.table("../../../Dataset/outs.csv", header = T, sep = ",")


params <- rbind(outputRegress[[1]][,54937:54999], outputRegress[[2]][,54937:54999], outputRegress[[3]][,54937:54999])
pred0 <- rbind(outputRegress[[1]][,1:18312], outputRegress[[2]][,1:18312], outputRegress[[3]][,1:18312])
pred1 <- rbind(outputRegress[[1]][,18313:36624], outputRegress[[2]][,18313:36624], outputRegress[[3]][,18313:36624])
pred2 <- rbind(outputRegress[[1]][,36625:54936], outputRegress[[2]][,36625:54936], outputRegress[[3]][,36625:54936])
Nin <- as.vector(as.matrix(ins))
Nout <- as.vector(as.matrix(outs))

#only the intercept
#Provare a sommare anche il termine della rain?

inter <- params[,seq(1,63, by=3)]
rain <- params[,seq(2,63, by=3)]
temp <- params[,seq(3,63, by=3)]

#We select a single lambda for example the correlation lambda_0
l0 <- c(1,2,7,8,9,10,19)
l1 <- c(3,4,11,12,13,14,20)
l2 <- c(5,6,15,16,17,18,21)

inter0 <- inter[,l0]
inter1 <- inter[,l1]

x11()
par(mfrow = c(4,4))
for(i in 1:4){
  for(j in 1:2){
    #Traceplot
    plot(inter[,7] + inter[,j] + inter[,2+i])
    #Density plot
    hist(inter[,7] + inter[,j] + inter[,2+i])
    
  }
}

#Check the predictions (maybe divide by chain(?))

Nin_pred <- pred0 + pred1
Nout_pred <- pred0 + pred2

alpha <- 0.025

Nin_q1 <- apply(Nin_pred, 2, function(x){quantile(x,alpha)})
Nin_q2 <- apply(Nin_pred, 2, function(x){quantile(x,1-alpha)})
Nout_q1 <- apply(Nout_pred, 2, function(x){quantile(x,alpha)})
Nout_q2 <- apply(Nout_pred, 2, function(x){quantile(x,1-alpha)})

covering_Nin <- sum(as.numeric(Nin_q1 <= Nin & Nin <= Nin_q2))/length(Nin)
covering_Nout <- sum(as.numeric(Nout_q1 <= Nout & Nout <= Nout_q2))/length(Nout)

covering_Nin 
covering_Nout 


#Plot a single node

nd <- 49

#Trace and acf plot of the node, first chain only


x11()
par(mfrow=c(1,2))
plot(exp(inter0[230:1000,1]+ inter0[230:1000,3]+inter0[230:1000,7])+exp(inter1[230:1000,1]+ inter1[230:1000,3]+inter1[230:1000,7]), type = "l", ylab = "Lambda of the N_in", main = "Traceplot")
acf(exp(inter0[230:1000,1]+ inter0[230:1000,3]+inter0[230:1000,7])+exp(inter1[230:1000,1]+ inter1[230:1000,3]+inter1[230:1000,7]), main = "Autocorrelation")


#Plot the confidence intervals for a single node over the time
#Changing the ylim might be necessary
pred_y <- Nin_pred[,((nd-1)*m + 1):(nd*m)]
pred <- colMeans(Nin_pred[,((nd-1)*m + 1):(nd*m)])

x11()
plot(1:m, pred, pch = 16, col = 'red', type ='o', xlab = 'day', ylab='traffic', main =paste("Prediction intervals for arriving bikes at cluster", nd, sep = " "), ylim = c(0,400))
points(1:m, Nin[((nd-1)*m + 1):(nd*m)], pch=15)
for(ii in 1:m)
{
  points(ii, quantile(pred_y[,ii], alpha), col = 'darkred', pch = 16)
  points(ii, quantile(pred_y[,ii], 1-alpha), col = 'darkred', pch = 16)
  segments(ii, quantile(pred_y[,ii] , alpha), ii, quantile(pred_y[,ii], 1-alpha), col  = 'darkred')
  if(ii!=42)
  {
    lines(c(ii, ii+1), c(quantile(pred_y[,ii] , alpha), quantile(pred_y[,ii+1], alpha)), col = 'darkred', lty = 3)
    lines(c(ii, ii+1), c(quantile(pred_y[,ii] , 1-alpha), quantile(pred_y[,ii+1], 1-alpha)), col = 'darkred', lty = 3)
  }
}
