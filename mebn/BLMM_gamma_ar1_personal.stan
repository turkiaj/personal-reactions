
data { 
  // new data
  int<lower=0> N;   // number of personal observations
  int<lower=1> p;   // number of predictors
  int<lower=1> k;   // number of group-level predictors
  matrix[N,p] X;    // fixed-effect design matrix
  matrix[N,k+1] Z;  // random-effect design matrix
  vector[N] Y;      // response

  // model parameters
  int<lower=0> N_samples;               // number of samples in parameter posteriors
  vector[N_samples] beta_Intercept;     // intercept 
  vector[p] beta[N_samples];            // poulation-level effects (fixed effects)
  vector[N_samples] ar1;                // AR(1) coefficient
  vector[N_samples] g_alpha;            // alpha (shape) parameter of the gamma distribution
  //matrix[k+1, k+1] Sigma_b;
  //matrix[k+1, k+1] L_prior;
  cholesky_factor_corr[k+1] L_prior[N_samples];
  vector<lower=0>[k+1] sigma_b_prior[N_samples];   // personal-effect standard deviations
} 

transformed data {
  vector[N] Y_t;      // transformed response
  
  for (n in 1:N)
  {
      Y_t[n] = Y[n] + 20;
  }  
}

parameters {
  //cholesky_factor_corr[k+1] L;    // Cholesky factor of group ranef corr matrix
  //vector<lower=0>[k+1] sigma_b;   // personal-effect standard deviations
  vector[k+1] z;                  // unscaled group-level effects
}

transformed parameters {
  vector[k+1] b;  

  for (s in 2000:N_samples)
  {
    b = diag_pre_multiply(sigma_b_prior[s], L_prior[s]) * z;   
  }
}

model { 
  real mu;              
  //vector[k+1] null_vector;
  
  // null_vector = rep_vector(0, k+1);
  z ~ normal(0,1);
  
  // sigma_b ~ student_t(3, 0, 10);
  // L ~ lkj_corr_cholesky(1);

  //target += multi_normal_cholesky_lpdf(b | null_vector, L_prior); 
  //target += multi_normal_lpdf(b | rep_vector(0, k+1), Sigma_b); 
  
  for (n in 1:N) 
  {
    for (s in 2000:N_samples)
    {
        mu = beta_Intercept[s] + 20 + X[n] * beta[s] + Z[n] * b;

        if (n > 1)
          mu += Y_t[n-1] * ar1[s];
  
        // Likelihood of the personal effect
        target += gamma_lpdf(Y_t[n] | g_alpha[s], g_alpha[s] / mu);
    }
  }
}

generated quantities { 
  vector[N] Y_rep;           // repeated response
  real mu_hat;

  //vector[k-1] personal_effect;

  // Posterior predictive distribution for model checking
  
  for (n in 1:N) 
  {
    for (s in 2000:N_samples)
    {
        mu_hat = beta_Intercept[s] + 20 + X[n] * beta[s] + Z[n] * b;

        if (n > 1)
          mu_hat += Y[n-1] * ar1[s];
  
        Y_rep[n] = gamma_rng(g_alpha[s], g_alpha[s] / mu_hat) - 20;
    }
  }

  // Personal effect combining typical effect and personal variation
  //personal_effect = beta + b[2:k+1];
}


