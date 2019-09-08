
data { 
  // new data
  int<lower=0> N;   // number of personal observations
  int<lower=1> p;   // number of fixed predictors
  int<lower=1> k;   // number of random predictors
  matrix[N,p] X;    // fixed-effect design matrix
  matrix[N,k] Z;    // random-effect design matrix
  vector[N] Y;      // response

  // model parameters
  int<lower=0> N_samples;               // number of samples in parameter posteriors
  matrix[N_samples, p] beta_Intercept;  // intercept 
  matrix[N_samples, p] beta;            // poulation-level effects (fixed effects)
  vector[N_samples] ar1;                // AR(1) coefficient
  vector[N_samples] g_alpha;            // alpha (shape) parameter of the gamma distribution
  matrix[k, k] Sigma_b[N_samples];  // variance-covariance matrix of group-level effects
} 

model { 
}

generated quantities { 
  matrix[N_samples, N] Y_pred;          // predicted response
  matrix[N_samples, p] personal_effect;
  vector[k] z;             // unscaled personal-level effects
  vector[k] b;               // personal-effects for predicted person
  real mu_hat;
  real g_beta_hat;

  for (n in 1:N)
  {
    for (s in 1:N_samples)
    {
      for (j in 1:k)
      {
        z[j] = normal_rng(0,1);
        
        b = Sigma_b[j] * z[j];
        personal_effect[s, j] = beta[s, j] + b[s, j];
      }
      
      mu_hat = beta_Intercept + X[n] * beta[s] + Z[n] * b[s];
      
      if (n > 1)
        mu_hat += Y[n-1] * ar1[s];
        
      g_beta_hat = g_alpha[s] / mu_hat;
      Y_pred[s, n] = gamma_rng(g_alpha[s], g_beta_hat);
      
    }
  }
} 
