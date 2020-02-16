#Plotting the Num(t) functions

rm(list=ls())
library(dplyr)
load("../../Dataset/bikemi_data.RData")
dat1<-read.table(file = "../../Dataset/hour.csv", sep = ",", header = T)
data <- cbind(bike_nil_week_hour, dat1)
rm(bike_nil_week_hour, dat1)

m <- 100  #number of intervals everyday
h <- 240000/m

t <- seq(0, 240000, by = h)

plot_flux <- function(station, day){
  sub <- data[which(data$NewDay == day), c(10,12, 24, 25)]
  sub_in <- sub[which(sub$NumeroStzFine==station), 4]
  sub_out <- sub[which(sub$NumeroStzInizio==station), 3]
  n_in <- rep(0, m+1)
  n_out <- rep(0, m+1)
  for(i in 2:(m+1)){
    n_in[i] <- as.numeric(length(which(sub_in < t[i])))
    n_out[i] <- as.numeric(length(which(sub_out < t[i])))
  }
  return(n_in - n_out)
}

plot(1:(m+1),plot_flux(2,3),xaxt='n',yaxt='n',xlab="time", ylab="Nin-Nout", 
     title(main="Net time-interval flow", sub="100 intervals in a day"), pch=19, cex=0.6,
     col="darkred")
axis(side=1,at=seq(0,100,by=4))
axis(side=2,at=seq(-100,100,by=10),las=1)
lines(1:(m+1),plot_flux(2,3),col="darkred")
