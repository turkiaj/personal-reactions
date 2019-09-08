
data { 
  int<lower=0> N;   // number of personal observations
  int<lower=0> N_samples;   // number of samples in parameter posteriors
  int<lower=1> p;   // number of predictors
  int<lower=1> k;   // number of group-level predictors
  matrix[N,p] X;    // fixed-effect design matrix
  matrix[N,k] Z;    // random-effect design matrix
  vector[N] Y;      // response

  real beta_Intercept;            // intercept 
  vector[p] beta;                 // poulation-level effects (fixed effects)
  real ar1;                       // AR(1) coefficient
  real<lower=0> g_alpha;          // alpha (shape) parameter of the gamma distribution
  matrix[k, k] Sigma_b;           // variance-covariance matrix of group-level effects
} 


model { 
  for (j in 1:k)
    z[j] ~ normal(0,1);
}

generated quantities { 
  vector[k] b;                   // personal-effects for predicted person
  vector[k] z;                   // unscaled group-level effects
  vector[N] Y_pred;               // predicted response

  b = Sigma_b * z;    

  group_size = 0;     // group variables for AR computation
  current_group = -1; // group variables for AR computation

  for (n in 1:N-NH)
  {
    mu_hat = beta_Intercept + X[n] * beta + Z[n] * b;
    
    if (current_group != group[n]) {
      current_group = group[n];
      group_size = 1;
    } 
    else
      group_size += 1;

    // - add autoregression coefficient from previous possibly correlated observation 
    if (group_size > 1)
      mu_hat += Y_t[n-1] * ar1;
      
      g_beta_hat = g_alpha / mu_hat;
      Y_rep[n] = gamma_rng(g_alpha, g_beta_hat);
  }
  
  // Personal effects for both training and holdout persons 
  for (j in 1:J) 
  {
    personal_effect[j] = beta + b[j][2:k];
  }  

} 
