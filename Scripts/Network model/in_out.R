#Node by node analysis

load("../../Dataset/bikemi_data.RData")
dat <- bike_nil_week_hour

stations <- levels(as.factor(dat$StazioneFine))
n <- length(stations)

ins <- as.data.frame(matrix(rep(0,42*n), 42, n))
names(ins) <- stations

outs <- as.data.frame(matrix(rep(0,42*n), 42, n))
names(outs) <- stations

for(i in 1:dim(dat)[1]){
  ins[dat[i, "NewDay"], dat[i, "StazioneFine"]] <- ins[dat[i, "NewDay"], dat[i, "StazioneFine"]] + 1
  outs[dat[i, "NewDay"], dat[i, "StazioneInizio"]] <- outs[dat[i, "NewDay"], dat[i, "StazioneInizio"]] + 1
  
}

write.csv(file = "../../Dataseet/ins.csv", ins, sep = ",", row.names = F)
write.csv(file = "../../Dataset/outs.csv", outs, sep = ",", row.names = F)

sapply(ins-outs, max)
sapply(ins-outs, min)