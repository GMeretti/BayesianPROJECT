#Prepares the datset for the network models
rm(list=ls())

#Decide how to divide the day
H <- 4  #number of time phases each day

load("../../Dataset/clusterized_data_M3.RData")
dat <- new_data_layer
rm(new_data_layer)

#hour <- read.table("../../Dataset/hour.csv", header = T, sep = ",")
time <- read.table("../../Dataset/time.csv", header = T, sep = ",")



stations <- levels(as.factor(dat$NumeroStzFine))
n <- length(stations)

#categorize by t,h over the time

#Uniform division over the day only when using the hour.csv dataset
dh <- (24*3600)/H

#time_ph <- hour%/%dh

time_ph <- time

ins <- as.data.frame(matrix(rep(0,42*n*H), 42*H, n))
names(ins) <- stations

outs <- as.data.frame(matrix(rep(0,42*n*H), 42*H, n))
names(outs) <- stations

#Loops over the data set and fills it up

for(i in 1:dim(dat)[1]){
  ins[as.numeric((dat[i, "NewDay"]-1)*H + time_ph[i,2]), as.character(dat[i, "NumeroStzFine"])] <- ins[as.numeric((dat[i, "NewDay"]-1)*H + time_ph[i,2]), as.character(dat[i, "NumeroStzFine"])] + 1
  outs[as.numeric((dat[i, "NewDay"]-1)*H + time_ph[i,1]), as.character(dat[i, "NumeroStzInizio"])] <- outs[as.numeric((dat[i, "NewDay"]-1)*H + time_ph[i,1]), as.character(dat[i, "NumeroStzInizio"])] + 1
  
}

write.csv(ins, file = "../../Dataset/ins.csv")
write.csv(outs, file = "../../Dataset/outs.csv")
