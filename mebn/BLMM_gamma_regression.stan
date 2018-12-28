
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
  vector[Pc] beta;                // population-level effects (fixed effects)
  real<lower=0> g_log_alpha;      // alpha (shape) parameter of the gamma distribution
}

transformed parameters {
  vector[N] mu;                  // local variable for mean of Normal distribution
  real<lower=0> g_alpha;         // alpha (shape) parameter of the gamma distribution
  vector<lower=0>[N] g_beta;     // beta (rate) of Gamma distribution
  real<lower=0> sigma_e;         // residual standard deviations 

  // - mean, or typical correlation
  mu = temp_Intercept + Xc * beta;
  
  // - log transform alpha parameter to keep it positive
  g_alpha = g_log_alpha;

  for (n in 1:N) 
  {
     // - use identity link for mu
     g_beta[n] = g_alpha / mu[n];
  }
  
  // variance can be estimated through transformation 
  // (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4024993/)
  sigma_e = log(1/g_alpha + 1);
}

model { 
  
  // Priors
  beta[1] ~ cauchy(0,10); // prior for the intercept following Gelman 2008

  for(i in 2:Pc)
   beta[i] ~ cauchy(0,2.5); //prior for the slopes following Gelman 2008
   
  // Likelihood with gamma regression
  Y ~ gamma(g_alpha, g_beta);
}

generated quantities { 
  real beta_Intercept;            // population-level intercept 
  vector[N] Y_rep;                // repeated response
  vector[N] log_lik;              // log-likelihood for LOO
  real mu_hat;
  real g_beta_hat;

  beta_Intercept = temp_Intercept - dot_product(means_X, beta);
  
  // Posterior predictive distribution for model checking

  for (n in 1:N) 
  {
    mu_hat = beta_Intercept + Xp[n] * beta;
    g_beta_hat = g_alpha / mu_hat;
    
    Y_rep[n] = gamma_rng(g_alpha, g_beta_hat);
    
    // Compute log-Likelihood for later LOO comparison of the models 
    log_lik[n] = gamma_lpdf(Y[n] | g_alpha, g_beta_hat);
  }
} 
