data{
  int<lower=0> N; //number of observations
  int J; //number of populations
  array[N] real spawners; 
  array[N] real survival; //log(Recruits/spawner)
  array[J] real Smax_mean; //prior mean of Smax
  array[J] real Smax_sigma; // prior sigma of Smax
  array[J] int start_row; //start of observation for that population
  array[J] int end_row; //end of observations for that population
  
}

transformed data {
  
  vector[J] log_Smax_mean;
  vector[J] log_Smax_sigma;
  
  //moment matching
  for (j in 1:J){
    log_Smax_mean[j] = log(Smax_mean[j]) - 0.5*log(1 + (Smax_sigma[j]^2)/(Smax_mean[j]^2)); //convert smax prior to per capita slope - transform to log scale with bias correction
    log_Smax_sigma[j] = sqrt(log(1 + (Smax_sigma[j]^2)/(Smax_mean[j]^2))); //this converts sigma on the untransformed scale to a log scale
  }
  
  
  
}

parameters {
  
  vector<lower=0>[J] Smax; // spawner level at which recruits are maximized (according to Ricker model)
  vector[J] alpha_j; //river level intrinsic productivity, intercept
  //variance components
  real<lower=0> mu_sigma; ///mean sigma among all stocks
}

transformed parameters {
  vector<lower=0>[J] b; //per capita density dependence term, slope
  //productivity residuals through time
  vector[N] e_t; //stock residual productivity at time t
  vector[N] mu; //expectation at each time for each stock
  
  for(j in 1:J){ // for every population
  b[j] = 1/Smax[j];
  
  for(t in start_row[j]:end_row[j]){ //for every year
  mu[t] = alpha_j[j] - b[j]*spawners[t];
  e_t[t] = survival[t] - mu[t]; // no autocorrelation in errors
  }
  
  }
  
}


model {
  alpha_j ~ normal(1.2,2); //prior for intrinsic productivity for all stocks
  
  for(j in 1:J){
    Smax[j] ~ lognormal(log_Smax_mean[j], log_Smax_sigma[j]); //prior on Smax for each stock
  }
  //variance terms
  mu_sigma ~ normal(1,1);
  
  //likelihood
  for(j in 1:J){
    survival ~ normal(mu, mu_sigma); //likelihood for all observations (no autocorrelation)
  }
  
  
}


generated quantities {
  
  
  
}
