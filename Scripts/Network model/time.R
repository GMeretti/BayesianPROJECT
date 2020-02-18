#Divides into the 4 time phases of the day, same used for the BSTS
#5-9 / 10-15 / 16-19/ 20-02

load("../../Dataset/bikemi_data.RData")
dat <- bike_nil_week_hour[, 22:23]
rm(bike_nil_week_hour)

limits <- c(9,15,19,24)

divide <- function(x){
  if(x == 1 | x == 2){return(4)}else{
    return(sum(as.numeric(x > limits)) + 1)
  }
}


time <- as.data.frame(matrix(0, nrow = 350093, ncol = 2))

time[] <- vapply(as.matrix(dat), divide, numeric(1))

write.csv(time, "../../Dataset/time.csv", row.names = F)
