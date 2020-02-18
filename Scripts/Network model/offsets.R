#Calculate the offsets 

#Zero information

write.csv(rep(1,109), file = "../../Dataset/offsets.csv", row.names = F)

#Taking information from a part of the dataset

ins <- read.table("../../Dataset/ins.csv", header = T, sep = ",")
outs <- read.table("../../Dataset/outs.csv", header = T, sep = ",")

#Use the first K time units to estimate the offsets

K <- 50

off <- (colMeans(ins[1:K,2:110]) + colMeans(outs[1:K,2:110]))/2

write.csv(off, file = "../../Dataset/offsets.csv", row.names = F)
