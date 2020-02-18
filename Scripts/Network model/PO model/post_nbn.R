#Analysis of the posterior

name <- "Poisson_nbn"

load(paste(name, "output.Rdata", sep = "_"))

ins <- read.table("../../../Dataset/ins.csv", header = T, sep = ",")
outs <- read.table("../../../Dataset/outs.csv", header = T, sep = ",")


params <- rbind(outputRegress[[1]][,c(54967:55002, 55112:55120)], outputRegress[[2]][,c(54967:55002, 55112:55120)], outputRegress[[3]][,c(54967:55002, 55112:55120)])
phi <- rbind(outputRegress[[1]][,55003:55111], outputRegress[[2]][,55003:55111], outputRegress[[3]][,55003:55111])
pred0 <- rbind(outputRegress[[1]][,37:18348], outputRegress[[2]][,37:18348], outputRegress[[3]][,37:18348])
pred1 <- rbind(outputRegress[[1]][,18349:36660], outputRegress[[2]][,18349:36660], outputRegress[[3]][,18349:36660])
pred2 <- rbind(outputRegress[[1]][,36661:54972], outputRegress[[2]][,36661:54972], outputRegress[[3]][,36661:54972])
Nin <- as.vector(as.matrix(ins))
Nout <- as.vector(as.matrix(outs))

#only the intercept

inter <- params[,seq(1,45, by=3)]
rain <- params[,seq(2,45, by=3)]
temp <- params[,seq(3,45, by=3)]

#We select a single lambda for example the correlation lambda_0
l0 <- c(13)
l1 <- c(1,2,5,6,7,8,14)
l2 <- c(3,4,9,10,11,12,15)

inter <- inter[l0]

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
m <- dim(ins)[1]

#Trace and acf plot of the node, first chain only
x11()
par(mfrow=c(1,2))
plot(exp(inter0[230:1000] + phi[230:1000,nd])+exp(inter1[230:1000,1]+ inter1[230:1000,3]+inter1[230:1000,7]+phi[230:1000,nd]), type = "l", ylab = "Lambda of the N_in", main = "Traceplot")
acf(exp(inter0[230:1000] + phi[230:1000,nd])+exp(inter1[230:1000,1]+ inter1[230:1000,3]+inter1[230:1000,7]+phi[230:1000,nd]), main = "Autocorrelation")

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
