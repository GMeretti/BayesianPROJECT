The files make use of the clusterized dataset M3

time.R : divides each observation in the 4 time phases
hour.R : converts time from strings to numeric for easier manipulation
offsets.R : create a dataset for the offsets used in the Offsets model

network_data.R : constructs the raw N_in N_out datasets 

analysis.R : preprocessing and checking the correlation

Model folders: each folder contains a .bug file with the model, a poi_*.R that builds it with JAGS and a post_*.R that analyse 
the posterior and predictive distributions

OLD
Poisson_try.bug
