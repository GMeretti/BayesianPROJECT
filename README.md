# Bayesian project

The repository is structured in the following way:
> Dataset: all dataset used for the analysis
--> bike_data.RData: original dataset with travels information
--> data_weather.csv: data concerning weather conditions, e.g. temperature (max,min,avg), umidity, rain, wind..
--> ins.csv / out.csv: matrix of inward / outward flows for the network model
--> hour.csv: exact times of departures and arrivals
--> locations_stations.csv: 
--> offset.csv: 
--> times.csv: 

> Presentations: the three presentation of the project

> Report: two reports are available, the short version (20 pages) and the complete version, to be taken as reference for any detail.

> Scripts: the scripts we implemented to perform the analysis. It is further subdivided into:
--> DBSCAN: related to the clustering processing phase, required for the network analysis
--> Frequentist stationarization: related to the stationarity analysis of the time series
--> Global Model: contains all the models concerning total volume analysis of the entire network
--> Network Model: contains all the scripts for the network approach
