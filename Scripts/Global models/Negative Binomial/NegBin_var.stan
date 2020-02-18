//#######################################################
//### GLMM - NegBin REGRESSION 	      #################
//#######################################################

///////////////////////// DATA /////////////////////////////////////
  data {
    int<lower = 0> N;       // number of data
    int<lower = 0> p_fix;   // number of covariates, fixed effects
    int<lower = 0> ngr;	// number of groups
    
    int<lower = 0> Y[N];  	// response vector
    matrix[N, p_fix] X;   	// design matrix (fixed effects)
    matrix[N, ngr] G;     	// groups (dummy) allocation
  }

//////////////////// PARAMETERS /////////////////////////////////
  parameters {
    vector[p_fix] beta;// regression coefficients 
    vector[ngr] r_param;
    vector[ngr] theta;      	// (group specific) random effects
    vector[p_fix] sigma2_beta;	// variances for the prior on beta
    vector[ngr] sigma2_theta;	// variances for the prior on theta
  }

//////////////////// TRANSFORMED PARAMETERS /////////////////////////////////
  transformed parameters 
{
  vector[N] mu;
  vector[N] r_val;
  for(i in 1:N){
    mu[i] = exp(row(X, i) * beta + row(G, i) * theta);
    r_val[i] = row(G, i) * r_param;
  }
}

////////////////// MODEL ////////////////////////
  model {
    
    // Likelihood     
    for (s in 1:N)
    {
      // print(" s = ", s);
      Y[s] ~ neg_binomial_2(mu[s], r_val[s]);  
    } 
    
    for (j in 1:p_fix) 
    {
      beta[j] ~ normal(0.0, pow(sigma2_beta[j], 0.5));
      sigma2_beta[j] ~ inv_gamma(2., 10.);
    }
    
    for (j in 1:ngr) 
    {
      theta[j] ~ normal(0.0, pow(sigma2_theta[j], 0.5));
      sigma2_theta[j] ~ inv_gamma(2., 10.);
    }
    
    for (j in 1:ngr) 
    {
      r_param[j] ~ uniform(0,50);
    }
    
  }

////////////////// GENERATED QUANTITIES ////////////////////////
  generated quantities 
{
  vector[N] log_lik;
  vector[N] Ypred;
  for (j in 1:N){
    log_lik[j] = poisson_lpmf(Y[j] | mu[j]);
    Ypred[j] ~ neg_binomial_2(mu[j], r_val[j]); 
  }
}
