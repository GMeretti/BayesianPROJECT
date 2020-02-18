rm(list=ls())

load("bikemi_data_torti.RData")
data <- bike_nil_week_hour             
head(data)

# day by day (subchoice)
attach(data)
data_wed <- data[which(data$DayOfTheWeek == 3),]
head(data_wed)

data_sun <- data[which(data$DayOfTheWeek == 7),]
head(data_sun)

data_sat <- data[which(data$DayOfTheWeek == 6),]
head(data_sat)

# cadorna -> duomo on wednesdays
CAD_DUO_wed <- data_wed[which(data_wed$StazioneInizio == "Cadorna 1" | data_wed$StazioneInizio == "Cadorna 2" | 
                                data_wed$StazioneInizio == "Cadorna 3"| data_wed$StazioneInizio == "Cadorna 4"),] 
CAD_DUO_wed <- CAD_DUO_wed[which(CAD_DUO_wed$StazioneFine == "Duomo"),]
head(CAD_DUO_wed)

unique(CAD_DUO_wed$Date1) # almost 13 travels a wednesday ( 78 / 6 wednesdays)
nrow(CAD_DUO_wed)

# cadorna -> duomo on sundays
CAD_DUO_sun <- data_sun[which(data_sun$StazioneInizio == "Cadorna 1" | data_sun$StazioneInizio == "Cadorna 2" | 
                                data_sun$StazioneInizio == "Cadorna 3"| data_sun$StazioneInizio == "Cadorna 4"),] 
CAD_DUO_sun <- CAD_DUO_sun[which(CAD_DUO_sun$StazioneFine == "Duomo"),]
head(CAD_DUO_sun)

unique(CAD_DUO_sun$Date1) # almost 1.3 travels a sunday ( 8 / 6 sundays (2 without any travel))
nrow(CAD_DUO_sun)

# cadorna -> duomo on saturdays
CAD_DUO_sat <- data_sat[which(data_sat$StazioneInizio == "Cadorna 1" | data_sat$StazioneInizio == "Cadorna 2" | 
                                data_sat$StazioneInizio == "Cadorna 3"| data_sat$StazioneInizio == "Cadorna 4"),] 
CAD_DUO_sat <- CAD_DUO_sat[which(CAD_DUO_sat$StazioneFine == "Duomo"),]
head(CAD_DUO_sat)

unique(CAD_DUO_sat$Date1) # almost 4.5 travels a satday ( 27 / 6 saturdays (1 without any travel))
nrow(CAD_DUO_sat)


# weekend vs weekdays

data_we <- data[which(data$week ==0),]
head(data_we)

data_wd <- data[which(data$week ==1),]
head(data_wd)

# cadorna -> duomo on weekends
CAD_DUO_we <- data_we[which(data_we$StazioneInizio == "Cadorna 1" | data_we$StazioneInizio == "Cadorna 2" | 
                                data_we$StazioneInizio == "Cadorna 3"| data_we$StazioneInizio == "Cadorna 4"),] 
CAD_DUO_we <- CAD_DUO_we[which(CAD_DUO_we$StazioneFine == "Duomo"),]
head(CAD_DUO_we)

nrow(CAD_DUO_we)
unique(CAD_DUO_we$Date1) # almost 3.1 (35/11 because the last sunday is incomplete) travels a sunday 
# ( 35 / 12 weekend days)

# cadorna -> duomo on weekdays
CAD_DUO_wd <- data_wd[which(data_wd$StazioneInizio == "Cadorna 1" | data_wd$StazioneInizio == "Cadorna 2" | 
                              data_wd$StazioneInizio == "Cadorna 3"| data_wd$StazioneInizio == "Cadorna 4"),] 
CAD_DUO_wd <- CAD_DUO_wd[which(CAD_DUO_wd$StazioneFine == "Duomo"),]
head(CAD_DUO_wd)

nrow(CAD_DUO_wd)
unique(CAD_DUO_wd$Date1) # almost 14.2 (426/30 because the last sunday is incomplete) travels a sunday 
# ( 426 / (5*6) weekdays )



# Duomo -> Cadorna on weekends
DUO_CAD_we <- data_we[which(data_we$StazioneFine == "Cadorna 1" | data_we$StazioneFine == "Cadorna 2" | 
                              data_we$StazioneFine == "Cadorna 3"| data_we$StazioneFine == "Cadorna 4"),] 
DUO_CAD_we <- DUO_CAD_we[which(DUO_CAD_we$StazioneInizio == "Duomo"),]
head(DUO_CAD_we)

nrow(DUO_CAD_we)
unique(DUO_CAD_we$Date1) # almost 2.4 (27/11 because the last sunday is incomplete) travels a sunday 
# ( 27/ 12 weekend days)

# Duomo -> Cadorna on weekdays
DUO_CAD_wd <- data_wd[which(data_wd$StazioneFine == "Cadorna 1" | data_wd$StazioneFine == "Cadorna 2" | 
                              data_wd$StazioneFine == "Cadorna 3"| data_wd$StazioneFine == "Cadorna 4"),] 
DUO_CAD_wd <- DUO_CAD_wd[which(DUO_CAD_wd$StazioneInizio == "Duomo"),]
head(DUO_CAD_wd)

nrow(DUO_CAD_wd)
unique(DUO_CAD_wd$Date1) # almost 13.7 (411/30 travels )
