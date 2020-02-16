##################################################
##### DBSCAN PREPROCESSING                   #####
##################################################
## Prerequisites
library(dbscan)
library(RColorBrewer)
library(dplyr)

## Load bike data
rm(list=ls())
load("../../Dataset/bikemi_data.RData")
data <- bike_nil_week_hour
rm(bike_nil_week_hour)

## Load locations data
locations <-read.table("../../Dataset/locations_stations.txt", header = T, sep = " ")
dim(locations)
N_loc <- dim(locations)[1] # total number of stations

####################################################
##### penalized DBSCAN method for clustering   #####
####################################################
## Definition of distances form the center function
#  --> @par locations: table of locations imported before
#  --> @pr x0: center coordinates
#  --> @ret a vector of distances from the selected center
distances_from_center <- function (pos, x0)
{
  n_dist <- dim(pos)[1]
  distances <- data.frame()
  for (i in 1:n_dist)
  {
    distances = rbind(distances, sqrt((pos$lat[i]-x0[1])^2+(pos$long[i]-x0[2])^2))
  }
  names(distances) = c("dist")
  return(distances)
}

## Penalization function
# --> @par gamma: in [0,1], limit for n->+infinity
# --> @par alpha: >= 0, so that penalization(gamma,1)==1,
#           it's a shape parameter, defautlt 1
# --> @par beta: in [0,1], trunction of the distance, default 0
# --> @ret the penalization evaluaed in the location
penalization <- function(x, alpha = 1, beta = 0, gamma)
{
  return(max(exp((x/alpha)^2*log(1-gamma))+gamma, beta))
}
##------------------------------------------------##

## Computing distances from the center and normalizing them
center =  c(45.458498166, 9.188165914)  # center are the exact locations of the dome
dist_ctr = distances_from_center(locations, center)  # vector of distances from the center
head(dist_ctr)

med = median(dist_ctr$dist) # mdian of the distances from the center

## REMARK
## From now on the median will represent the scale factor for the
## distances in the Milan network: all distances will be normalized via this term
## for a more intuitive interpretation of the penalization parameters

summary_distances = c(min = min(dist_ctr$dist),    qt05 = quantile(dist_ctr$dist, 0.05),
                      med = median(dist_ctr$dist), qt95 = quantile(dist_ctr$dist, 0.95),
                      max = max(dist_ctr$dist))
summary_distances

x11()
hist(dist_ctr$dist, breaks = 20)  # histogram of the distances

# Normalization of diatnces using the median
normalized_distances = data.frame(dist = dist_ctr$dist/median(dist_ctr$dist))
normalized_distances

summary_ndistances = c(min = min(normalized_distances$dist),    qt05 = quantile(normalized_distances$dist, 0.05),
                       med = median(normalized_distances$dist), qt95 = quantile(normalized_distances$dist, 0.95),
                       max = max(normalized_distances$dist))
summary_ndistances

x11()
hist(normalized_distances$dist, breaks = 20) # plot of the normalized distances distribution
graphics.off()                               # it's obviously just a rescaling of the previous one
##------------------------------------------------##

## Definition of function semidistance
# --> @par gamma, alpha: are the shape parameters of the penalization
# --> @par x0: is the city center
# --> @par locations: locations dataset
# --> @par nor_dist: table of normalized distances frm the center
# --> @par med: median of the original distribution of distances from the center
# --> @ret semidistance matrix of the dataset
semidistance <- function(alpha, beta, gamma, x0, pos, nor_dist, med, with_idx = FALSE, idx)
{
  n_matrix <- dim(pos)[1]
  semi_matrix <- matrix(rep(0,n_matrix*n_matrix), n_matrix, n_matrix)
  for(i in 1:n_matrix)
  {
    if(with_idx == FALSE)
    {
      dis_i = nor_dist[i,1]
    } else {
      dis_i = nor_dist[idx[i],1]
    }
    pen_i = penalization(dis_i, alpha, beta, gamma)
    for(j in 1:n_matrix)
    {
      if(with_idx == FALSE)
      {
        dis_j = nor_dist[j,1]
      } else {
        dis_j = nor_dist[idx[j],1]
      }
      pen_j = penalization(dis_j, alpha, beta, gamma)
      pen = (pen_i+pen_j)/2
      dis = sqrt((pos$lat[i]-pos$lat[j])^2+(pos$long[i]-pos$long[j])^2)/med
      semi_matrix[i,j] = dis*pen
    }
  }
  return(semi_matrix)
}
##------------------------------------------------##

## LOOP loss function
# --> @par locations: locations dataset
# --> @par cluster: clustering partition of the nodes
# --> @ret percentage of oops in the cluster
loop_loss <- function(trips, pos, cldiv, sumq)
{
  ## Data augmentation (try it with new and old model)
  #  -> dataset + last column dbscan value
  scanned_stations = data.frame(pos, db = cldiv+sumq)
  indices = which(scanned_stations$db==sumq)
  scanned_stations[indices, ]$db = scanned_stations[indices, ]$id
  head(scanned_stations)
  
  augmented_data = data.frame(start = trips$NumeroStzInizio, end = trips$NumeroStzFine)
  augmented_data = augmented_data %>% left_join(scanned_stations, by=c("start"="id")) %>% rename(db_start = db)
  augmented_data = data.frame(db_start=augmented_data$db_start, end=augmented_data$end)
  augmented_data = augmented_data %>% left_join(scanned_stations, by=c("end"="id")) %>% rename(db_end = db)
  augmented_data = data.frame(db_start=augmented_data$db_start, db_end=augmented_data$db_end)
  
  ## Percentage of induced loops
  ## Augment the corresponding selected model to see results
  loops <- augmented_data %>% filter(db_start==db_end)
  n_loops = dim(loops)[1]
  
  perc_loops = n_loops/dim(data)[1]
  return(perc_loops)
}

## LOOP PERCENTAGE IN VARIOUS GROUPS FUNCTION
print_loop_distribution <- function (trips, cutoff)
{
  loops_grouped <- trips %>% filter(NumeroStzInizio==NumeroStzFine) %>%
    group_by(NumeroStzInizio) %>% summarize(n())
  grouped <- trips %>% group_by(NumeroStzInizio) %>% summarize(n())
  perc = loops_grouped[,2] /grouped[,2]
  
  # Big custers
  print("Big clusters:")
  big_cl = which(grouped[,2]>10000)
  data_big = data.frame(grouped[big_cl,], percentage=perc[big_cl,])
  print(data_big)
  
  # Outliers
  print("Idices outliers:")
  outl = which(perc[,1]>cutoff)
  print(outl)
  
  print("Outliers:")
  data_out = data.frame(grouped[outl,], percentage=perc[outl,])
  print(data_out)
  
  # Histogram with outlier
  x11()
  hist(perc[,1], breaks = 20, col='lightblue', 
       main = 'Histogram of loop percentages', xlab='Loop percentage')
  # Histogram without outliers
  x11()
  hist(perc[-outl,1], breaks = 18, col='lightblue', 
       main = 'Histogram of loop percentages', xlab='Loop percentage')
  # Histogram of numbers
  x11()
  hist(as.data.frame(grouped[,2])[,1], breaks=20, col='lightblue',
       xlab='Number of out-links', main='Distribution of out-links')
  
  print("IC:")
  IC = cbind(q05 = quantile(perc[,1], 0.05), med = median(perc[,1]),
             qt95 = quantile(perc[,1],0.95))
  print(IC)
  
  print("IC without outliers:")
  IC_without = cbind(q05 = quantile(perc[-outl,1], 0.05), med = median(perc[-outl,1]),
             qt95 = quantile(perc[-outl,1],0.95))
  print(IC_without)
}
##------------------------------------------------##

## Parameters selection in DBSCAN algorithm
## -> ep is th epsilon distance
## -> alpha, gamma used in the penalization
evaluation <- data.frame()
for(alpha in seq(from=0.60, by=0.05, to=0.90))
{
  for(gamma in seq(from=0.5, by=0.05, to=0.6))
  {
    for (ep in seq(from=0.90, by=0.05, to=1.15))
    {
      mtx <- semidistance(alpha, 0, gamma, center, locations, normalized_distances, med)
      cluster <- dbscan(mtx, eps = ep, minPts = 2)  ## produces clusters and noise-points
      print(paste0("CLUSTER(ep = " , ep, ", alpha = ", alpha,
                   ", gamma = ", gamma, "): ", length(unique(cluster$cluster))))   ## number_of_clusers+1 (noise points are labeled as 0)
      leng <- length(unique(cluster$cluster))
      sing_perc <- table(cluster$cluster)[["0"]]/N_loc
      loop_perc <- loop_loss(data, locations, cluster$cluster, 1000)
      evaluation <- rbind(evaluation, cbind(alpha, gamma, ep, leng, sing_perc, loop_perc))
      
      #if(leng>30) ## not to few clusters
      #{
      #  x11()
      #  plot(locations[,2], locations[,1], pch = ((cluster$cluster %% 25)+(cluster$cluster>=25)*(cluster$cluster!=0)),
      #       xlab = 'longitude', ylab = 'latitude',
      #       main=paste('epsilon = ',ep, ", alpha = ", alpha, ", gamma = ", gamma),
      #       col=rainbow(100)[locations$nil+10], ylim=c(45.44 , 45.537), xlim = c(9.08, 9.3))
      #  legend(9.25, 45.53, legend = unique(locations$nilname), cex = 0.5,
      #         col=rainbow(100)[unique(locations$nil)+10], pch = 16)
      # Colors are Milan neighbourhoods, symbols are clusters, empty square means singleton
      #}
    }
  }
}
names(evaluation) = c("alpha", "gamma", "epsilon", "number", "singletons_perc", "loop_perc")
evaluation

## Discard cases where minimal gof assumptions are not met
evaluation_refined <- evaluation %>% filter(loop_perc<0.15, singletons_perc<0.33, number>40)
evaluation_refined

## Find the best model with an automatic evaluation
selected_model = evaluation_refined[which.min(evaluation_refined$singletons_perc+3*evaluation_refined$loop_perc), ]
selected_model

## Plot selected model
mtx = semidistance(0.8, 0, 0.5, center_duomo, locations, normalized_distances, med)
cluster <- dbscan(mtx, eps = 1.15, minPts = 2)

x11()
plot(locations[,2], locations[,1], pch = cluster$cluster, xlab = 'longitude', ylab = 'latitude',
     main=paste('Proposed clustering'),
     col=rainbow(100)[locations$nil+10], ylim=c(45.44 , 45.537), xlim = c(9.08, 9.3))
legend(9.25, 45.53, legend = unique(locations$nilname), cex = 0.5,
       col=rainbow(100)[unique(locations$nil)+10], pch = 16)
dev.off()
##------------------------------------------------##

## Disaggregating big clusters
disag = 13

dt = locations[which(cluster$cluster == disag),]
head(dt)

alpha = 0.8
beta  = 0
gamma = 0.5

idx = which(cluster$cluster == disag)
semi_matrix <- semidistance(alpha, beta, gamma, x0, dt, normalized_distances, med, TRUE, idx)

evaluation <- data.frame()
for (ep in seq(from=0.1, by=0.05, to=1.15))
{
  cl <- dbscan(semi_matrix, eps = ep, minPts = 2)  ## produces clusters and noise-points
  print(paste0("CLUSTER(ep = " , ep, ", alpha = ", alpha,
               ", gamma = ", gamma, "): ", length(unique(cl$cluster))))   ## number_of_clusers+1 (noise points are labeled as 0)
  leng <- length(unique(cl$cluster))+ (sum(cl$cluster == 0) -1*as.numeric(sum(cl$cluster == 0)!=0))
  perc_loops = loop_loss(data, dt, cl$cluster, 2000)
  
  evaluation <- rbind(evaluation, cbind(alpha, gamma, ep, leng, perc_loops))
  
  x11()
  plot(dt[,2], dt[,1], pch = ((cl$cluster %% 25)+(cl$cluster>=25)*(cl$cluster!=0)),
       xlab = 'longitude', ylab = 'latitude',
       main=paste('epsilon = ',ep, ", alpha = ", alpha, ", gamma = ", gamma),
       col=rainbow(100)[dt$nil+10], ylim=c(45.44 , 45.537), xlim = c(9.08, 9.3))
  legend(9.25, 45.53, legend = unique(locations$nilname), cex = 0.5,
         col=rainbow(100)[unique(locations$nil)+10], pch = 16)
}
graphics.off()

names(evaluation) = c("alpha", "gamma", "epsilon", "number", "loop_perc")
evaluation

## best
ep = 0.35
cl <- dbscan(semi_matrix, eps = ep, minPts = 2)

normalized_loop_percentage = evaluation$loop_perc/max(evaluation$loop_perc)
normalized_loop_percentage

x11()
plot(seq(from=0.1, by=0.05, to=1.15),normalized_loop_percentage, type='b', xlab='eps', ylab='Normalized loop percentage in cluster 1013', ylim=c(0,1))
x11()
plot(seq(from=0.1, by=0.05, to=1.15), evaluation$number, type='b', xlab='eps', ylab='Number of subgroups in cluster 1013')
graphics.off()
##------------------------------------------------##

## OLD MODEL OF LAST YEAR
## Parameters selection in DBSCAN algorithm
## -> ep is th epsilon distance
## -> mp is the minumum number of points in a cluster
for (ep in seq(from=0.001, by=0.001, to=0.01))
{
  for(mp in 2:5)
  {
    cluster <- dbscan(locations[ ,1:2], eps = ep, minPts = mp)  ## produces clusters and noise-points
    print(paste0("CLUSTER(ep = " , ep, ", mp = ", mp, "): ", length(unique(cluster$cluster))))   ## number_of_clusers+1 (noise points are labeled as 0)
    if(length(unique(cluster$cluster))>10) ## not to few clusters
    {
      x11()
      plot(locations[,2], locations[,1], pch = cluster$cluster, xlab = 'longitude', ylab = 'latitude',
           main=paste('epsilon = ',ep, ' MinPoints = ', mp), 
           col=rainbow(100)[locations$nil+10], ylim=c(45.44 , 45.537), xlim = c(9.08, 9.3))
      legend(9.25, 45.53, legend = unique(locations$nilname), cex = 0.5,  
             col=rainbow(100)[unique(locations$nil)+10], pch = 16)
      ## Colors are Milan neighbourhoods, symbols are clusters, empty square means singleton
    }
  }
}
graphics.off()

## Model selected last year
dist_table = dist(locations[,1:2])
cluster <- dbscan(dist_table, eps = 0.004, minPts = 2)

x11()
plot(locations[,2], locations[,1], pch = cluster$cluster, xlab = 'longitude', ylab = 'latitude',
     main=paste('epsilon = ', 0.04, ' MinPoints = ', 2), 
     col=rainbow(100)[locations$nil+10], ylim=c(45.44 , 45.537), xlim = c(9.08, 9.3))
legend(9.25, 45.53, legend = unique(locations$nilname), cex = 0.5,  
       col=rainbow(100)[unique(locations$nil)+10], pch = 16)

table(cluster$cluster)  ## high disomogeneity
dev.off()
##------------------------------------------------##

#### LOOP EVALUATION
## Origianl amount of loops in dataset
or_loops <- data %>% filter(NumeroStzInizio==NumeroStzFine)
or_n_loops = dim(or_loops)[1]

or_perc_loops = or_n_loops/dim(data)[1]
or_perc_loops

print_loop_distribution(data, 0.06)
graphics.off()
##------------------------------------------------##

####################################################
#####   dataset creation with selected model   #####
####################################################
## REMARK: Before doing this remember to have as
## cluster the mdoel you want to save (!!!)
cldiv = cluster$cluster

## !!! Correction if needed !!!
cldiv[idx[which(cl$cluster==2)]] = length(table(cluster$cluster)) # add new subcluster

cldiv

##--------------------------------------------##

#### USEFUL FUNCTIONS
build_new_stations <- function(pos, cldiv, sumq)
{
  scanned_stations = data.frame(pos, db = cldiv)
  head(scanned_stations)
  
  ## calculation of new cluster baricentrs
  baricenters <- scanned_stations %>% filter(db>0) %>% group_by(db) %>% summarise(lat = mean(lat), long = mean(long)) %>% ungroup()
  baricenters$db <- baricenters$db+sumq ## creating a new id different from the old ones
  baricenters = data.frame(baricenters[,c(2,3,1)])
  colnames(baricenters) = c('lat', 'long', 'id')
  zeros <- data.frame(lat = scanned_stations$lat[which(scanned_stations$db==0)],
                      long = scanned_stations$long[which(scanned_stations$db==0)],
                      id = scanned_stations$id[which(scanned_stations$db == 0)]) ## filtering the singletons
  new_stations <- rbind(zeros, baricenters) ## adding the new clusters to the filtered singletons

  return(new_stations)
}
##------------------------------------------------##

## Computing NETWORK FLOW
## Definition of conversion matrix: translates the old id in
## the new cluster based one, if there are singletons
## the id is unchanged
build_new_dataset <- function(trips, pos, cldiv, sumq, layer = FALSE)
{
  N_loc = dim(pos)[1]
  conversion <- data.frame(id = pos$id, db = cldiv)
  conversion$db <- conversion$db+sumq
  
  for(i in 1:N_loc)
  {
    if(conversion$db[i]==sumq)
      conversion$db[i] = conversion$id[i]
  }
  
  M <- dim(trips)[1]
  
  new_data <- trips
  if(layer == FALSE)
  {
    new_data <- new_data %>% left_join(conversion, by=c("NumeroStzInizio"="id")) %>%
      rename(OldNumeroStzInizio = NumeroStzInizio) %>% rename(NumeroStzInizio = db) %>%
      left_join(conversion, by=c("NumeroStzFine"="id")) %>% rename(OldNumeroStzFine = NumeroStzFine) %>%
      rename(NumeroStzFine = db)
  } else {
    new_data <- new_data %>% left_join(conversion, by=c("NumeroStzInizio"="id")) %>%
      rename(DoubleOldNumeroStzInizio = OldNumeroStzInizio) %>%
      rename(OldNumeroStzInizio = NumeroStzInizio) %>% rename(NumeroStzInizio = db) %>%
      left_join(conversion, by=c("NumeroStzFine"="id")) %>% rename(DoubleOldNumeroStzFine = OldNumeroStzFine) %>%
      rename(OldNumeroStzFine = NumeroStzFine) %>% rename(NumeroStzFine = db)
  }
  
  return(new_data)
}
##------------------------------------------------##

## crete the new dataset
new_data = build_new_dataset(data, locations, cldiv, 1000)
head(new_data)

new_stations = build_new_stations(locations, cldiv, 1000)
head(new_stations)

N_stat <- dim(new_stations)[1] ## number of selected post-dbscan clusters (macro-stations)
##------------------------------------------------##

## To save the file
addendum = "M2"
save (new_stations, file = paste0("../../Dataset/clusterized_locations_", addendum, ".RData"))
save (new_data, file = paste0("../../Dataset/clusterized_data_", addendum, ".RData"))
##------------------------------------------------##

#### EVALUATION
## LOOP LOSS
## new/old model, use old data to force augmentation
perc_loops <- loop_loss(data, locations, cldiv, 1000)
perc_loops

## Loop percentage distribution
print_loop_distribution(new_data, 0.10)
graphics.off()
##------------------------------------------------##

####################################################
#####         second layer clustering          #####
####################################################
evaluation_l <- data.frame()

alpha = 0.8
gammma = 0.5

for (ep in seq(from=0.90, by=0.05, to=1.15))
{
  mtx <- semidistance(alpha, 0, gamma, center, new_stations, normalized_distances, med)
  cluster <- dbscan(mtx, eps = ep, minPts = 2)  ## produces clusters and noise-points
  print(paste0("CLUSTER(ep = " , ep, ", alpha = ", alpha,
               ", gamma = ", gamma, "): ", length(unique(cluster$cluster))))   ## number_of_clusers+1 (noise points are labeled as 0)
  leng <- length(unique(cluster$cluster))
  sing_perc <- table(cluster$cluster)[["0"]]/N_loc
  loop_perc <- loop_loss(data, new_stations, cluster$cluster, 2000)
  evaluation_l <- rbind(evaluation_l, cbind(alpha, gamma, ep, leng, sing_perc, loop_perc))
  
  #if(leng>30) ## not to few clusters
  {
    x11()
    plot(new_stations[,2], new_stations[,1], pch = ((cluster$cluster %% 25)+(cluster$cluster>=25)*(cluster$cluster!=0)),
        xlab = 'longitude', ylab = 'latitude',
         main=paste('epsilon = ',ep, ", alpha = ", alpha, ", gamma = ", gamma), ylim=c(45.44 , 45.537), xlim = c(9.08, 9.3))
  }
}
graphics.off()

names(evaluation_l) = c("alpha", "gamma", "epsilon", "number", "singletons_perc", "loop_perc")
evaluation_l

## Discard cases where minimal gof assumptions are not met
evaluation_refined_l <- evaluation_l %>% filter(loop_perc<0.15)
evaluation_refined_l

## Find the best model with an automatic evaluation
selected_model_l = evaluation_refined_l[which.min(evaluation_refined_l$singletons_perc+3*evaluation_refined_l$loop_perc), ]
selected_model_l

## best
ep_l = 1.1

mtx_layer <- semidistance(alpha, 0, gamma, center, new_stations, normalized_distances, med)
cl_layer <- dbscan(mtx_layer, eps = ep_l, minPts = 2)

x11()
plot(new_stations[,2], new_stations[,1], pch = cl_layer$cluster, xlab = 'longitude', ylab = 'latitude',
     main=paste('Proposed clustering'), ylim=c(45.44 , 45.537), xlim = c(9.08, 9.3))
graphics.off()
##------------------------------------------------##

## crete the new dataset
new_data_layer = build_new_dataset(new_data, new_stations, cl_layer$cluster, 2000, TRUE)
head(new_data_layer)

new_stations_layer = build_new_stations(new_stations, cl_layer$cluster, 2000)
head(new_stations_layer)
##------------------------------------------------##

## To save the file
addendum = "M3"
save (new_stations_layer, file = paste0("../../Dataset/clusterized_locations_", addendum, ".RData"))
save (new_data_layer, file = paste0("../../Dataset/clusterized_data_", addendum, ".RData"))
##------------------------------------------------##

#### EVALUATION
## LOOP LOSS
## new/old model, use old data to force augmentation
perc_loops_layer <- loop_loss(new_data, new_stations, cl_layer$cluster, 2000)
perc_loops_layer

## Loop percentage distribution
print_loop_distribution(new_data_layer, 0.10)
graphics.off()

## Monitoring clusters
num = 1
first_level = new_stations[which(cl_layer$cluster == num),]$id
up_lv = first_level[which(first_level>1000)]-1000
down_lv = first_level[which(first_level<1000)]
locations[c(which(cldiv == up_lv), down_lv),]

####################################################
#####      some useful extra information       #####
####################################################
## Definition of flow matrix
N_stat = ...
adj_mtx <- matrix(rep(0,N_stat^2), N_stat, N_stat)
db_names <- new_stations$db

for (i in 1:dim(new_data)[1])
{
  idi <- db_names==new_data$Start_db[i]
  idj <- db_names==new_data$End_db[i]
  adj_mtx[idi,idj] <- adj_mtx[idi,idj] + 1
}

for (i in 1:67)
  adj_mtx[i,i] <- 0 ## auto-connections are not considered

rownames(adj_mtx)<-db_names
colnames(adj_mtx)<-db_names

adj_mtx

write.table(adj_mtx, "../../Dataset/Extra/flow_matrix.txt")
##------------------------------------------------##

## Computation of GEOODESIC DISTANCE matrix
library(geosphere)

geo_dist <-matrix(rep(0, N_stat*N_stat), N_stat)

for(i in 1:N_stat)
{
  for(j in 1:N_stat)
  { # distances comuted in meters, the matrix is symmetric
    geo_dist[i,j]<-distGeo(c(new_stations[i,2], new_stations[i,3]),
                           c(new_stations[j,2], new_stations[j,3])) 
    
  }
}

colnames(geo_dist) <- new_stations$db
rownames(geo_dist) <- new_stations$db
head(geo_dist)

write.table(geo_dist, file = "../../Dataset/Extra/geodesic_distances.txt")
