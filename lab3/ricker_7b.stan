data{
  int<lower=0> N; //number of observations
  int J; //number of populations
  int K; // number of exposure categories
  int A; // number of areas
  int years; // numer of years
  array[N] int ii; //index of brood years
  array[N] int aa; //index of areas
  array[N] real spawners; 
  array[N] real survival; //log(Recruits/spawner)
  matrix[N,K] exposure_category_matrix; // dummy variable with 1,0
  array[J] int start_row; //start of observation for that population
  array[J] int end_row; //end of observations for that population
  
}



parameters {
  
  
  vector[years] theta_year; // random effect of year
  real<lower=0> sigma_theta_area; // standard deviation for the distribution of random effects theta_year_area 
  vector[K] beta; //effect of exposure on survival
  
  //variance components
  real<lower=0> mu_sigma; // standard deviation for survival
  vector[years] z; // parameter that affects theta_year_area (also called non centered parameterization)
  vector<lower=0>[J] b; //per capita density dependence term, slope
}

transformed parameters {
  
  vector[N] exposure; //total exposure effect for each observation
  vector[N] mu; //expectation at each time for each stock
  matrix[A,years] theta_year_area; //nested random effect of area within year
  
  for(a in 1:A){ # loop over areas
    for(year in 1:years){ # loop over years
      theta_year_area[a,year]  = theta_year[year] +  z[year]*sigma_theta_area; //nested random effect of area within year
      //distribution of theta_year_area effect has a mean which is theta_year
    }
    
  }
  
  
  for(j in 1:J){ // for every population
  
  for(t in start_row[j]:end_row[j]){ //for every year
  exposure[t] = 0; // set exposure to 0
  for(k in 1:K){// for every level in the exposure category
    exposure[t] = exposure[t] + exposure_category_matrix[t,k]*beta[k]; 
    // sum over the exposures. exposure category matrix has 1,0. beta is the effect of exposure on survival
    
  }
  
  mu[t] = exposure[t] + theta_year_area[aa[t],ii[t]]  - b[j]*spawners[t]; // survival = r + theta -bS
  }
  
  }
  
}


model {
  
  theta_year ~ normal(0,1); // prior for random effect of year  
  b ~ normal(0,0.1); // prior for slope
  //variance terms
  mu_sigma ~ normal(1,1); // prior for standard deviation
  sigma_theta_area ~ normal(1,1); // prior for standard deviation for random effect
  
  //likelihood
  for(j in 1:J){ //for every river
    survival ~ normal(mu, mu_sigma); //likelihood for all observations 
  }
  
  
}




