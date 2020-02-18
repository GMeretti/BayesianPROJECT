#Hour time modifications
#Used when dividing the day in equal intervals
rm(list=ls())
load("../../Dataset/bikemi_data.RData")
data <- bike_nil_week_hour
rm(bike_nil_week_hour)

hour1 <- rep(0, dim(data)[1])
hour2 <- rep(0, dim(data)[1])
hour <- as.data.frame(cbind(hour1, hour2))

#Transforms the strings in numeric format, number of seconds from midnight

transform_hour <- function(st){
  M <- strsplit(st, " ")[[1]][2]
  aux <- strsplit(st, ":")[[1]]
  aux[3] <- strsplit(aux[3], " ")[[1]][1]
  res <- sum(as.numeric(aux)*c(3600, 60, 1))
  if(((M == "p.") & (res%/%3600 != 12)) | ((M == "a.") & (res%/%3600 == 12))){
    res <- res + 12*3600
  }
  return(res)
}

for(i in 1:dim(data)[1]){
  hour[i,1] <- transform_hour(as.character(data$Hour1[i]))
  hour[i,2] <- transform_hour(as.character(data$Hour2[i]))
  
}

write.csv(hour, file = "../../Dataset/hour.csv", row.names = F)
