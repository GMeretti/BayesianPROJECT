#Analysis of ins and outs

ins <- read.table("../../Dataset/ins.csv", header = T, sep = ",")
outs <- read.table("../../Dataset/outs.csv", header = T, sep = ",")

#We divide the ins and outs inside the 4 time phases

cv <- matrix(0, nrow = 109, ncol = 4)
cr <- matrix(0, nrow = 109, ncol = 4)

for(k in 1:4){
  for(i in 1:109){
    cv[i,k] <- cov(ins[seq(k,168, by=4),i], outs[seq(k,168, by=4),i])
    cr[i,k] <- cor(ins[seq(k,168, by=4),i], outs[seq(k,168, by=4),i])
  }
}

#Eliminate NAs

cr <- cr[-47,]

#Correlation plot

x11()

par(mfrow=c(2,2))
for(i in 1:4){
  plot(cr[,i], pch=19, xlab = "Nodes", ylab = "Correlation", main = paste0("Phase", i))
  abline(h=mean(cr[,i]), col="red", lwd=4)
}


#Plot of a single node

x11()
nd <- 73
week <- rep(c(rep(2,20),rep(4,8)), each=6)

par(mfrow=c(1,2))
plot(ins[,nd],outs[,nd], col = rep(1:4, 42), main = "Divided by time phase", xlab = "N in", ylab = "N out", pch=19)
plot(ins[,nd],outs[,nd], col = week, main = "Divided by weekend/weekday" , xlab = "N in", ylab = "N out", pch=19)



