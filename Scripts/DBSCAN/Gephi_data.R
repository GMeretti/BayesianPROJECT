##################################################
##### GEPHI                                  #####
##################################################

###### Original dataset
## Prerquisites
library(dplyr,tidyr)
library(rgexf)

rm(list=ls())
load("../../Dataset/bikemi_data.RData")
data <- bike_nil_week_hour

locations <-read.table("../../Dataset/locations_stations.txt", header = T, sep = " ")
dim(locations)
N_loc <- dim(locations)[1] # total number of stations

#### Gephi representation adapted form already existing material
## Define GLOBAL network
network <- data %>% group_by(NumeroStzInizio,NumeroStzFine) %>% summarise(Flow = n()) %>% 
  ungroup() %>% rename(Source = NumeroStzInizio, Target = NumeroStzFine)
head(network)

## Define GLOBAL nodes
nodes <- locations %>% select(id, nome, lat, long)
head(nodes)

## Modify names to make them readable from Gephi
geonodes <- nodes%>%select(id,nome)
geonodes$nome <- as.character(geonodes$nome)

geonodes$nome[47]  <- "Santissima Trinita"  # corrections for safe encoding
geonodes$nome[91]  <- "San Lorenzo - Colonne"
geonodes$nome[151] <- "Cantu"
geonodes$nome[180] <- "Donizetti Provincia - RIMOSSA PROVVISORIAMENTE"
geonodes$nome[198] <- "Universita Bocconi"
geonodes$nome[217] <- "M. Gioia - Regione"
geonodes$nome[235] <- "Universita Cattolica"

## Select longitude and latitude as attributes
geonodes_attr <- nodes %>% select(lat,long)
names(geonodes_attr) <- c("Latitude", "Longitude")

## Define edges of the graph [[Note very clear what is thE difference with original....]]
edges <- network %>% left_join(geonodes, by=c("Source"="id")) %>% select(Source,Target,Flow) %>% rename(IDS=Source)
edges <- edges   %>% left_join(geonodes, by=c("Target"="id")) %>% select(IDS,Target,Flow)    %>% rename(IDT=Target)
head(edges)

## WRITE FILE
write.gexf(nodes = geonodes, edges = edges[ ,1:2], edgesWeight = edges$Flow,
           nodesAtt = geonodes_attr, defaultedgetype="directed", output="../../Dataset/Gephi/original_net.gexf")
##------------------------------------------------##

###### Clusterized dataset
## Prerquisites
library(dplyr,tidyr)
library(rgexf)

rm(list=ls())
load("../../Dataset/clusterized_data_M2.RData")
data <- new_data

load("../../Dataset/clusterized_locations_M2.RData")
locations <- new_stations
dim(locations)
N_loc <- dim(locations)[1] # total number of stations

#### Gephi representation adapted form already existing material
## Define GLOBAL network
network <- data %>% group_by(NumeroStzInizio,NumeroStzFine) %>% summarise(Flow = n()) %>% 
  ungroup() %>% rename(Source = NumeroStzInizio, Target = NumeroStzFine)
head(network)

## Define GLOBAL nodes
nodes <- locations %>% select(id = id, lat = lat, long = long)
nodes_name = as.character(nodes$id)
head(nodes)

## Modify names to make them readable from Gephi
geonodes <- nodes%>%select(id)
geonodes = cbind(geonodes, nodes_name)
head(geonodes)

## Select longitude and latitude as attributes
geonodes_attr <- nodes %>% select(lat,long)
names(geonodes_attr) <- c("Latitude", "Longitude")

## Define edges of the graph [[Note very clear what is thE difference with original....]]
edges <- network %>% left_join(geonodes, by=c("Source"="id")) %>% select(Source,Target,Flow) %>% rename(IDS=Source)
edges <- edges   %>% left_join(geonodes, by=c("Target"="id")) %>% select(IDS,Target,Flow)    %>% rename(IDT=Target)
head(edges)

## WRITE FILE
write.gexf(nodes = geonodes, edges = edges[ ,1:2], edgesWeight = edges$Flow,
           nodesAtt = geonodes_attr, defaultedgetype="directed", output="../../Dataset/Gephi/M2_net.gexf")
