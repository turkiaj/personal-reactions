
data { 
  int<lower=0> N;   // number of observations
  int<lower=1> p;   // number of predictors
  int<lower=1> J;   // number of groups in data (persons)
  int<lower=1> k;   // number of group-level predictors
  int<lower=1,upper=J> group[N]; //group indicator
  matrix[N,p] X;    // fixed-effect design matrix
  matrix[N,k] Z;    // random-effect design matrix
  vector[N] Y;      // response
  real offset;
  
  // Holdout data

  int<lower=0> N_new;       // number of holdout observations
  matrix[N_new,p-1] X_new;  // fixed-effect design matrix (no intercept)
  matrix[N_new,k] Z_new;    // random-effect design matrix (no intercept)
  vector[N_new] Y_new;      // response
  int<lower=1> J_new;       // number of groups in holdout data (persons)
  int<lower=1,upper=J_new> group_new[N_new]; //group indicator
  int<lower=0,upper=1> predict_with_holdout; // are these values used for prediction?
  
  // Horseshoe prior data
  real<lower=0> scale_icept;    // prior std for the intercept
  real<lower=0> scale_global;   // scale for the half -t prior for tau
  real<lower=1> nu_global;      // degrees of freedom for the half -t priors for tau
  real<lower=1> nu_local;       // degrees of freedom for the half - t priors for lambdas
  real<lower=0> slab_scale;     // slab scale for the regularized horseshoe
  real<lower=0> slab_df;        // slab degrees of freedom for the regularized horseshoe
} 

transformed data { 
  // Centering data for more stable sampling 
  int Pc; 
  matrix[N, p - 1] Xc;    // X centered
  matrix[N, p - 1] Xp;    // X without intercept, non-centered
  vector[p - 1] means_X;  // column means of X before centering 
  
  Pc = p - 1;  // the intercept is removed from the design matrix 
  for (i in 2:p) { 
     means_X[i - 1] = mean(X[, i]); 
     Xc[, i - 1] = X[, i] - means_X[i - 1]; 
     Xp[, i - 1] = X[, i];
  }
}

parameters { 
  real temp_Intercept;            // temporary intercept 
  cholesky_factor_corr[k] L;      // Cholesky factor of group ranef corr matrix
  vector<lower=0>[k] sigma_b;     // group-level random-effect standard deviations
  vector[k] z[J];                 // unscaled group-level effects
  vector[k] z_new[J_new];         // unscaled group-level effects
  real<lower=0> g_log_alpha;      // alpha (shape) parameter of the gamma distribution

  real ar1;                       // AR(1) coefficient

  // horseshoe shrinkage parameters 
  real <lower=0> tau;             // global shrinkage parameter
  vector <lower=0>[Pc] lambda;    // local shrinkage parameter
  real<lower=0> caux;
  vector[Pc] zbeta;    
}

transformed parameters {
  real<lower=0> g_alpha;          // alpha (shape) parameter of the gamma distribution
  vector<lower=0>[N] g_beta;      // beta (rate) of Gamma distribution
  matrix[k, k] Sigma_b;           // variance-covariance matrix of group-level effects
  real<lower=0> sigma_e;          // residual standard deviations 
  vector[Pc] beta;                // poulation-level effects (fixed effects)
  vector[k] b[J];                 // group-level effects (random effects)
  
  vector[N] e;                    // residuals

  // Regularized horseshoe parameters
  real<lower=0> c;                    // slab scale 
  vector<lower=0>[Pc] lambda_tilde;   // 'truncated' local shrinkage parameter 

  // Apply regularizing horseshoe prior for betas
  c = slab_scale * sqrt(caux);
  lambda_tilde = sqrt(c^2 * square(lambda) ./ (c^2 + tau^2* square(lambda)));
  beta = zbeta .* lambda_tilde * tau;
  
  // Premultiply diagonal matrix [sigma_b] with the Cholesky decomposition L of
  // the correlation matrix Sigma_b to get variance-covariance matrix of group-level effects

  // diag(sigma_b) * L
  Sigma_b = diag_pre_multiply(sigma_b, L); 
  
  // Group-level effects are generated by multipying D (Sigma_b) with z 
  // that has standard normal distribution
    
  for(j in 1:J) 
    b[j] = Sigma_b * z[j];    
    
  // - log transform alpha parameter to keep it positive
  g_alpha = g_log_alpha;

  {  // scope for local variables
    vector[N] mu; // local variable for mean of Normal distribution

    // group variables for AR computation
    int group_size = 0;
    int current_group = -1;

    // - mean, or typical correlation
    mu = offset + temp_Intercept + Xc * beta;
    
    for (n in 1:N) 
    {
       if (current_group != group[n]) {
         current_group = group[n];
         group_size = 1;
       } 
       else
         group_size += 1;
  
       // - add group effects
       mu[n] += Z[n] * b[group[n]]; 
  
      // - residuals     
      e[n] = Y[n] - mu[n]; 
  
      // - add autoregression coefficient from previous observation 
      if (group_size > 1)
        mu[n] += e[n-1] * ar1;

      g_beta[n] = g_alpha / mu[n];
    }  
  } // local scope
  
  // - use identity link for mu
  //g_beta = g_alpha / mu;

  // estimate of variance 
  // (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4024993/)
  sigma_e = log(1/g_alpha + 1);
}

model { 
  
  // Finnish Horseshoe (Regularized Horseshoe) prior
  // half-t priors for lambdas and tau, and inverse-gamma for c^2
  
  zbeta ~ normal(0, 1); 
  lambda ~ student_t(nu_local, 0, 1);
  tau ~ student_t(nu_global, 0, scale_global * sigma_e);
  caux ~ inv_gamma(0.5 * slab_df, 0.5 * slab_df);
  
  sigma_b ~ student_t(3, 0, 10);
  L ~ lkj_corr_cholesky(1);
  
  // Standard normal prior for random effects
  for (j in 1:J)
    z[j] ~ normal(0,1);

  // Predicted random effects
  for (j in 1:J_new)
    z_new[j] ~ normal(0,1);

  // Likelihood: hierarchical gamma regression 

  Y ~ gamma(g_alpha, g_beta);
}

generated quantities { 
  real beta_Intercept;            // population-level intercept 
  corr_matrix[k] C;               // correlation matrix 
  vector[N] Y_rep;                // repeated response
  vector[N] log_lik;              // log-likelihood for LOO
  real mu_hat;
  real g_beta_hat;

  vector[N_new] Y_pred;           // repeated response
  vector[k] b_pred[J_new];        // group-level effects (random effects) for new patients
  vector[k-1] personal_effect[J_new];

  // Correlation matrix of random-effects, C = L'L
  C = multiply_lower_tri_self_transpose(L); 
  
  //beta_Intercept = temp_Intercept - dot_product(means_X, beta) - offset;
 beta_Intercept = temp_Intercept + offset;
 
  // Posterior predictive distribution for model checking

  for (n in 1:N) 
  {
      mu_hat = beta_Intercept + Xp[n] * beta + Z[n] * b[group[n]];
      g_beta_hat = g_alpha / mu_hat;
      Y_rep[n] = gamma_rng(g_alpha, g_beta_hat);
      
      // Compute log-Likelihood for later LOO comparison of the models 
      log_lik[n] = gamma_lpdf(Y[n] | g_alpha, g_beta_hat);
  }
  
  if (predict_with_holdout == 1) {
    
    // Personal effects for new patients are generated by multipying D (Sigma_b) with z 
    // that has standard normal distribution
    
    for(j in 1:J_new) 
      b_pred[j] = Sigma_b * z_new[j];    

    for (n in 1:N_new) 
    {
        mu_hat = beta_Intercept + X_new[n] * beta + Z_new[n] * b_pred[group_new[n]];
        g_beta_hat = g_alpha / mu_hat;
        Y_pred[n] = gamma_rng(g_alpha, g_beta_hat);
    }
    
    // sample personal effects for new patients
    for (j in 1:J_new) 
    {
      personal_effect[j] = beta + b_pred[j][2:k];
    }
  }  
} 
