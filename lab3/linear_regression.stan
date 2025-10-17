
data {
  int<lower=0> N; //number of observations or years in this case
  vector[N] S; //predictor of size N
  vector[N] y; //response variable of size N
}
parameters {
  real r; //intercept to be estimated
  real b; //slope to be estimated
  real<lower=0> sigma; //variance to be estimated
}
model {

  r ~ normal(1.2, 2); //prior for intercept, based on literature
  b ~ normal(0, 10); //prior for slope
  sigma ~ normal(0, 1); //prior for variance
  y ~ normal(r + b * S, sigma); // linear model
}

