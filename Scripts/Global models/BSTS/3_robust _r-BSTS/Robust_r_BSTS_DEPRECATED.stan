data
{
	int<lower = 0> sz;
	int<lower = 0> S;
	vector<lower = 0, upper = 20000>[sz] Y;

	real mean_mu;
	vector[S-1] mean_gamma;
	matrix[S-1,S-1] var_gamma;

	int R[sz];
	int T[sz];
}

parameters
{
	real avg_mu[sz];
	real avg_delta[sz];
	real rho[sz];
	real mu[sz];
	real delta[sz];
	real gamma[sz+S-2];
	real alpha;
	real beta[2];
	real sigma[5];
}

model
{
	// Prior:
	for (k in 1:5)
	{
		sigma[k] ~ uniform(50,400);
	}
	for (h in 1:2)
	{
		beta[h] ~ normal(0, 1000);
	}
	mu[1] ~ normal(mean_mu, 2*100);
	delta[1] ~ normal(0, 2*100);
	rho[1] ~ normal(0, 2*100);
	for (j in 1:(S-1))
	{
		gamma[j] ~ normal(mean_gamma[j], var_gamma[j,j]);
	}

	// Likelihood:
	for(i in 1:(sz-1))
	{
		if(i>2*S) //fake derterminism, low variance
		{
			avg_mu[i] ~ normal(0.5 *mu[i] + 0.33*mu[i-S] + 0.17*mu[i-2*S], 0.01);
			avg_delta[i] ~ normal(0.34*delta[i] + 0.33*delta[i-S] + 0.33*delta[i-2*S], 0.01);
		} else if(i>S) {
			avg_mu[i] ~ normal(0.66*mu[i] + 0.34*mu[i-S], 0.01);
			avg_delta[i] ~ normal(0.5*delta[i] + 0.5*delta[i-S], 0.01);
		} else {
			avg_mu[i] ~  normal(mu[i], 0.01);
			avg_delta[i] ~  normal(delta[i], 0.01);
		}

		mu[i+1] ~ normal(avg_mu[i] + avg_delta[i], 10000);
		delta[i+1] ~ normal(avg_delta[i], 10000);
		gamma[S-1+i] ~ normal(-sum(gamma[i:(i+S-2)]), 10000);
		rho[i+1] ~ normal(alpha*rho[i], 10000);
		Y[i] ~ normal(mu[i] + gamma[i+S-2] + rho[i] + beta[1]*R[i] + beta[2]*T[i], 10000);
	}

	Y[sz] ~ normal(mu[sz] + gamma[sz+S-2] + rho[sz] + beta[1]*R[sz] + beta[2]*T[sz], 10000);
}
