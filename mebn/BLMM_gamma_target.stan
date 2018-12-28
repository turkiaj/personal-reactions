
data { 
  int<lower=0> N;   // number of observations
  int<lower=1> p;   // number of predictors
  int<lower=1> J;   // number of groups in data (persons)
  int<lower=1> k;   // number of group-level predictors
  int<lower=1,upper=J> group[N]; //group indicator
  matrix[N,p] X;    // fixed-effect design matrix
  matrix[N,k] Z;    // random-effect design matrix
  vector[N] Y;      // response
}

parameters { 
  vector<lower=0>[N] g_alpha; // alpha (shape) parameter of the gamma distribution
  vector<lower=0>[N] g_beta;  // beta (rate) of Gamma distribution
}

model {
    Y ~ gamma(g_alpha, g_beta);
}

generated quantities { 
  vector[N] Y_rep;                // repeated response

  // Posterior predictive distribution for model checking

  for (n in 1:N) 
  {
    Y_rep[n] = gamma_rng(g_alpha[n], g_beta[n]);
  }
} 
