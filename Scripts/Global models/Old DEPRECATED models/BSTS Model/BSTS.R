#################################################
#####              BSTS models              #####
#################################################
## Prerequisites
rm(list=ls())
library(dplyr)
library(bsts)
load("../../../Dataset/bikemi_data.RData")
data <- bike_nil_week_hour
head(data)

weather_data <- read.table("../../../Dataset/data_weather.csv", header = T, sep = ",")
head(weather_data)

total_traffic <- data %>% group_by(NewDay) %>% summarize(n())
names(total_traffic) <- c("Day","N_travels")
##---------------------------------------------##

#### BSTS models
## Model 1: state space is only trend(mu_t) + seasonal(gamma_t)
# Define state space:
ss <- AddLocalLinearTrend(list(), total_traffic$N_travels)
ss <- AddSeasonal(ss, total_traffic$N_travels, nseasons = 7)
model1 <- bsts(total_traffic$N_travels,
               state.specification = ss,
               niter = 10000)

x11()
plot(model1)
x11()
plot(model1, "components")
graphics.off()

## Model 2: state space is only trend(mu_t) + seasonal(gamma_t) + complete_regression
# Define full dataset:
complete_data = cbind(total_traffic, weather_data$Tmean, weather_data$rain, weather_data$Wind_mean)
complete_data = complete_data[,2:5]
head(complete_data)

model2 <- bsts(N_travels ~ ., state.specification = ss, niter = 100000, data = complete_data)

x11()
plot(model2)
x11()
plot(model2, "comp")
x11()
plot(model2, "coef") # inclusion probabilities
graphis.off()

## Model 3: state space is only trend(mu_t) + seasonal(gamma_t) + partial_regression
# Define partial dataset:
partial_data = cbind(total_traffic, tmean = weather_data$Tmean, rain=weather_data$rain)
partial_data = partial_data[,2:4]
head(partial_data)

model3 <- bsts(N_travels ~ ., state.specification = ss, niter = 100000, data = partial_data)

x11()
plot(model3)
x11()
plot(model3, "comp")
x11()
plot(model3, "coef")
graphics.off()

## Compare the first 3 models
CompareBstsModels(list("Model 1" = model1,
                       "Model 2" = model2,
                       "Model 3" = model3),
                 colors = c("black", "red", "blue"))
dev.off()

## Try simple predictions of a week with model1
pred1 <- predict(model1, horizon = 7)

x11()
plot(pred1, plot.original = 42)
dev.off()

## Model 4: state space is only trend(mu_t) + seasonal(gamma_t) + partial_regression + Ar(1)
# Define new ss:
ss1 = AddAr(ss, lags = 2, y = partial_data$N_travels)
model4 <- bsts(N_travels ~ ., ss1, data = partial_data, niter = 20000)

x11()
plot(model4)
x11()
plot(model4, "comp")
x11()
plot(model4, "coef")
graphics.off()

## Model 5: model 4 - temperature
model5 <- bsts(N_travels ~ rain, state.specification = ss, niter = 100000, data = partial_data)

x11()
plot(model5)
x11()
plot(model5, "comp")
x11()
plot(model5, "coef")
graphics.off()

## Final comparison
CompareBstsModels(list("Model 1" = model1,
                       "Model 2" = model2,
                       "Model 3" = model3,
                       "Model 4" = model4,
                       "Model 5" = model5),
                 colors = c("black", "red", "blue", "green", 'orange'))
 dev.off()