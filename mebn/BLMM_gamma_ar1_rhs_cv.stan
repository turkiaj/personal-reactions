
data { 
  int<lower=0> N;   // number of all observations
  int<lower=0> NH;  // number of holdout observations
  int<lower=1> p;   // number of predictors
  int<lower=1> J;   // number of groups in data (persons)
  int<lower=1> k;   // number of group-level predictors
  int<lower=1,upper=J> group[N]; //group indicator
  matrix[N,p] X;    // fixed-effect design matrix
  matrix[N,k] Z;    // random-effect design matrix
  vector[N] Y;      // response
  real offset;      // constant temporarily added to intercept to keeps values positive 
  
  // indicates if observation is used for training (0) or prediction (1)
  int<lower=0,upper=1> holdout[N]; 
  
  // Parameters for Regularized Horseshoe Prior
  real<lower=0> scale_icept;    // prior std for the intercept
  real<lower=0> scale_global;   // scale for the half -t prior for tau
  real<lower=1> nu_global;      // degrees of freedom for the half -t priors for tau
  real<lower=1> nu_local;       // degrees of freedom for the half - t priors for lambdas
  real<lower=0> slab_scale;     // slab scale for the regularized horseshoe
  real<lower=0> slab_df;        // slab degrees of freedom for the regularized horseshoe
} 

transformed data { 

  matrix[N-NH,p-1] X_t; // training input
  matrix[N-NH,k] Z_t;   // training input
  vector[N-NH] Y_t;     // training response
  int t=1;              // index 

  matrix[NH,p-1] X_h;   // holdout input
  matrix[NH,k] Z_h;     // holdout input
  vector[NH] Y_h;       // holdout response
  int h=1;              // index 
  
  for (n in 1:N)
  {
    if (holdout[n] == 1)
    {
      // the intercept is removed from the design matrix 
      X_h[h,1:p-1] = X[n,2:p];
      Z_h[h] = Z[n];
      Y_h[h] = Y[n] + offset;
      h += 1;
    }
    else
    {
      // the intercept is removed from the design matrix 
      X_t[t,1:p-1] = X[n,2:p];
      Z_t[t] = Z[n];
      Y_t[t] = Y[n] + offset;
      t += 1;
    }
  }  
}

parameters { 
  real beta_Intercept;            // intercept 
  cholesky_factor_corr[k] L;      // Cholesky factor of group ranef corr matrix
  vector<lower=0>[k] sigma_b;     // group-level random-effect standard deviations
  vector[k] z[J];                 // unscaled group-level effects
  real<lower=0> g_log_alpha;      // alpha (shape) parameter of the gamma distribution
  real ar1;                       // AR(1) coefficient

  // horseshoe shrinkage parameters 
  real <lower=0> tau;             // global shrinkage parameter
  vector <lower=0>[p-1] lambda;   // local shrinkage parameter
  real<lower=0> caux;
  vector[p-1] zbeta;    
}

transformed parameters {
  real<lower=0> g_alpha;          // alpha (shape) parameter of the gamma distribution
  real<lower=0> sigma_e;          // residual standard deviations 
  vector[p-1] beta;               // poulation-level effects (fixed effects)
  vector[k] b[J];                 // group-level effects (random effects)
  
  // Regularized horseshoe parameters
  real<lower=0> c;                    // slab scale 
  vector<lower=0>[p-1] lambda_tilde;  // 'truncated' local shrinkage parameter 
  matrix[k,k] Lambda;

  // Apply regularizing horseshoe prior for betas
  c = slab_scale * sqrt(caux);
  lambda_tilde = sqrt(c^2 * square(lambda) ./ (c^2 + tau^2* square(lambda)));
  beta = zbeta .* lambda_tilde * tau;

  Lambda = diag_pre_multiply(sigma_b, L); 

  // Group-level effects are generated by multipying D (Sigma_b) with z 
  // that has standard normal distribution
    
  for(j in 1:J) 
    b[j] = Lambda * z[j];    
    
  // - log transform alpha parameter to keep it positive
  g_alpha = g_log_alpha;

  // estimate of variance 
  // (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4024993/)
  sigma_e = log(1/g_alpha + 1);
}

model { 
  real mu;                // linear prediction 
  real g_beta;            // beta (rate) of Gamma distribution
  int group_size = 0;     // group variables for AR computation
  int current_group = -1; // group variables for AR computation
  
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
  
  // Gelman 2008
  beta_Intercept ~ cauchy(0,10);

  // Calculate posterior probability with training data only
  for (n in 1:N-NH) 
  {
    // typical effect
    mu = beta_Intercept + offset + X_t[n] * beta;
    
    // - add personal effects
    mu += Z_t[n] * b[group[n]];
  
    if (current_group != group[n]) {
      current_group = group[n];
      group_size = 1;
    } 
    else
      group_size += 1;
  
    // - add autoregression coefficient from previous possibly correlated observation
    if (group_size > 1)
      mu += Y_t[n-1] * ar1;

    // - use identity link for mu
    g_beta = g_alpha / mu;

    // Likelihood 
    target += gamma_lpdf(Y_t[n] | g_alpha, g_beta);
  }
}

generated quantities { 
  vector[NH] Y_pred;              // predicted response
  real mu_hat;
  real g_beta_hat;
  int group_size = 0;     // group variables for AR computation
  int current_group = -1; // group variables for AR computation

  vector[k-1] personal_effect[J];

  // Posterior prediction for holdout person
  
  for (n in 1:NH)
  {
    // - offset constant is needed to model transformed data correctly
    mu_hat = beta_Intercept + offset + X_h[n] * beta + Z_h[n] * b[group[n]];
    
    if (current_group != group[n]) {
      current_group = group[n];
      group_size = 1;
    } 
    else
      group_size += 1;

    // - add autoregression coefficient from previous possibly correlated observation 
    if (group_size > 1)
      mu_hat += Y_h[n-1] * ar1;
      
      g_beta_hat = g_alpha / mu_hat;
      
      // - the offset constant is substracted from final prediction
      Y_pred[n] = gamma_rng(g_alpha, g_beta_hat) - offset;
  }

  // Personal effects for both training and holdout persons 
  for (j in 1:J) 
  {
    personal_effect[j] = beta + b[j][2:k];
  }
} 