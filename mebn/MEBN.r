#
# Functions for constructing
# Mixed Effects Bayesian Network
#
# Jari Turkia
#

##################################################

mebn.new_graph_with_randomvariables <- function(datadesc)
{
  library(igraph)  

  # Initialize graph structure
  reaction_graph <- make_empty_graph() 
  
  # TODO: Get nodes styles from description, and not like this
  min_order <- min(datadesc$Order)
  predictor_columns <- datadesc[datadesc$Order==min_order,]
  
  next_order <- max(datadesc$Order)
  assumed_targets <- datadesc[datadesc$Order==next_order,]
  
  # TODO: Get shape?
  
  # Add nodes to reaction graph for all the random variables  
  reaction_graph <- reaction_graph + vertices(as.vector(assumed_targets$Name), 
                                              label=as.vector(assumed_targets$Description), 
                                              type=next_order,
                                              color = "#74aaf2", 
                                              size = 1,
                                              shape = "circle")
                                              
                                              # shape = "distbox",
                                              # mean=50,
                                              # l95CI=30,
                                              # u95CI=70,
                                              # scalemin=0,
                                              # scalemax=100)
  
  reaction_graph <- reaction_graph + vertices(as.vector(predictor_columns$Name), 
                                              label=as.vector(predictor_columns$Description), 
                                              type=min_order,
                                              color = "#3cd164", 
                                              size = 1, 
                                              shape = "circle")
  
                                              # shape = "distbox",
                                              # mean=50,
                                              # l95CI=30,
                                              # u95CI=70,
                                              # scalemin=0,
                                              # scalemax=100)
  

  return(reaction_graph)
}

##################################################

mebn.posterior <- function(igraph_node)
{

}

##################################################
  
mebn.add_priornodes <- function(datadesc, reaction_graph)
{
  require(igraph)  
  
  # Add prior information from datadesc for random variables
  # - filter prior information for the nutrition predictors
  
  predictor_columns <- datadesc[datadesc$Order==100,]
  
  predictor_names <- as.vector(predictor_columns$Name)
  predprior <- datadesc[(!is.na(datadesc$Lowerbound) | !is.na(datadesc$Upperbound)) & datadesc$Order == 100,]
  predprior$PriorName <- paste0(predprior$Name, "_prior")
  
  # - add prior nodes
  reaction_graph <- reaction_graph + vertices(as.vector(predprior$PriorName), 
                                              label=as.vector(predprior$PriorName), 
                                              color = "#92d3a5", 
                                              size = 1, 
                                              shape = "distbox",
                                              mean=as.vector(predprior$Lowerbound+(predprior$Upperbound-predprior$Lowerbound)/2),
                                              l95CI=as.vector(predprior$Lowerbound),
                                              u95CI=as.vector(predprior$Upperbound),
                                              scalemin=as.vector(predprior$ScaleMin),
                                              scalemax=as.vector(predprior$ScaleMax))
  
  # Connect priors to predictors
  for (p in 1:nrow(predprior))
  {
    print(paste0(predprior[p,]$PriorName, " -> ", predprior[p,]$Name))
    reaction_graph <- reaction_graph + edges(as.vector(predprior[p,]$PriorName), as.vector(predprior[p,]$Name), shape = "arrow", weight = 1)
  }
  
  return(reaction_graph)
}

##################################################

mebn.normalize <- function(x, org_min, org_max) { return ((x - org_min) / (org_max - org_min)) }

##################################################

mebn.renormalize <- function(x, org_min, org_max) { return ((x + org_min) * (org_max - org_min)) }

##################################################

mebn.scale <- function(x, org_sd) 
{ 
  #n_scaled <- ((x - org_mean) / org_sd) + 100
  n_scaled <- x / org_sd
  return(n_scaled)
}

##################################################

mebn.rescale <- function(x, org_sd) 
{ 
  #n_rescaled <- ((x * org_sd) + org_mean)
  n_rescaled <- x * org_sd
  return(n_rescaled)
}

##################################################

mebn.scale_gaussians <- function(r, data, datadesc) 
{ 
  # data parameters contains whole dataset 
  # - subset only field(s) in datadesc so that index r matches
  
  data <- data[as.vector(datadesc$Name)]
  s <- data[,r] # data from column r
  
  if (datadesc[r,]$Distribution == "Gaussian")
  {
    # TODO: Scale also returns the scaling factor. Store it to restore the original scale.
    s <- mebn.scale(s, sd(s))
    #s <- mebn.normalize(s, min(s), max(s))
  }
  
  return (s) 
}

##################################################

mebn.set_model_parameters <- function(predictor_columns, target_column, group_column, inputdata, normalize_values, reg_params = NULL)
{
  ident <- function(x) { return (x) }
  predictors <- inputdata[as.vector(predictor_columns$Name)]
  
  target_name <- as.vector(target_column$Name)
  
  prior_sigma <- rep(-1, nrow(predictor_columns))
  dim(prior_sigma) <- nrow(predictor_columns)
  prior_mean <- rep(0, nrow(predictor_columns))
  dim(prior_mean) <- nrow(predictor_columns)
  
  # Set informative priors 
  if (!is.na(predictor_columns$Lowerbound) && !is.na(predictor_columns$Upperbound))
  {  
    prior_sigma <- c(predictor_columns$Upperbound - predictor_columns$Lowerbound)
    prior_mean <- c(predictor_columns$Lowerbound + prior_sigma / 2)
    
    dim(prior_sigma) <- length(prior_sigma)
    dim(prior_mean) <- length(prior_mean)
  }
  
  # Scale if the predictor is Gaussian
  N <- nrow(inputdata)
  
  if (normalize_values == TRUE)
  {
    X <- sapply(1:nrow(assumedpredictors), mebn.scale_gaussians, data = inputdata, datadesc = assumedpredictors)
    
    # append intercept 
    X <- cbind(rep(1,N), X)
    
    #Y <- scale(inputdata[target_name][,], center = FALSE, scale = TRUE)[,1]
    #Y <- mebn.scale(inputdata[target_name], sd(inputdata[target_name][,]))[,]  
    #Y <- mebn.normalize(inputdata[target_name], min(inputdata[target_name][,]), max(inputdata[target_name][,]))[,]  
    Y <- inputdata[target_name][,]
  }
  else
  {
    X <- cbind(rep(1,N), apply(predictors, 2, ident))
    Y <- inputdata[target_name][,]
  }

  params <- within(list(),
                   {
                     N <- N
                     X <- X
                     p <- k <- ncol(X)               # all predictors may have random effects
                     # Mean and variance of Gaussian prior predictors
                     X_prior_sigma <- prior_sigma    # no prior for the intercept 
                     X_prior_mean <- prior_mean      # sigma < 0 means noninformative prior
                     Y <- Y
                     Z <- X                      
                     J <- length(levels(inputdata[[group_column]]))
                     group <- as.integer(inputdata[[group_column]])
                     offset <- 15
                   })
  
  params <- c(params, reg_params)
  
  return(params)
}

##################################################

mebn.localsummary <- function(fit)
{
  draws <- extract(fit)
  
  # TODO: get CI from stan
  #fixef_CI <- apply(draws$beta, 2, quantile, probs = c(0.05, 0.95)) # 2 = by column

  ModelSummary <- within(list(),
                         {
                           intmean     <- mean(draws$beta_Intercept)
                           fixef       <- colMeans(draws$beta)
                           fixef_l95CI <- 0
                           fixef_u95CI <- 0
                           ranef_sd    <- colMeans(draws$sigma_b)
                           std_error   <- mean(draws$sigma_e)
                           Sigma_b     <- colMeans(draws$Sigma_b)
                           C           <- colMeans(draws$C)
                         })
  
  # Create matrix D
  sdM <- diag(ModelSummary$ranef_sd)
  ModelSummary$D <- sdM %*% ModelSummary$C %*% t(sdM)
  
  return(ModelSummary)
}

##################################################

mebn.get_localfit <- function(target_name, local_model_cache = "models")
{
  modelcache <- paste0(local_model_cache, "/", target_name, "_blmm", ".rds")
  print(paste0("Loading ", modelcache))
  
  localfit <- NULL
  
  if (file.exists(modelcache))
  {
    localfit <- readRDS(modelcache)
  }
  
  return(localfit)
}

##################################################

mebn.get_parents <- function(g, nodename)
{
  vindex <- as.numeric(V(g)[nodename])
  parents <- neighbors(g, vindex, mode = c("in"))
  
  return(parents)
}

##################################################

mebn.get_parents_with_type <- function(g, nodename, type)
{
  vindex <- as.numeric(V(g)[nodename])
  parents <- neighbors(g, vindex, mode = c("in"))
  
  # filter by type
  parents <- parents[parents$type==type]
  
  return(parents)
}

##################################################

mebn.LOO_comparison <- function(target_variables, graphdir1, graphdir2)
{
  library(loo)
  
  comparison<-data.frame(matrix(nrow=nrow(target_variables), ncol=3))
  colnames(comparison) <- c("distribution", graphdir1, graphdir2)
  
  n <- 1
  for (targetname in target_variables$Name)
  {
    # Get models to compare
    m1 <- mebn.get_localfit(paste0(graphdir1, "/", targetname))
    m2 <- mebn.get_localfit(paste0(graphdir2, "/", targetname))
    
    # Statistics for model 1
    m1_loglik <- extract_log_lik(m1, merge_chains = FALSE)
    
    remove(m1_rel_n_eff)
    remove(m1_loo)
    
    if (!any(is.na(exp(m1_loglik))))
    {
      m1_rel_n_eff <- relative_eff(exp(m1_loglik))
      suppressWarnings(m1_loo <- loo(m1_loglik, r_eff = m1_rel_n_eff, cores = 4))
    }
    
    # Statistics for model 2
    m2_loglik <- extract_log_lik(m2, merge_chains = FALSE)
    
    remove(m2_rel_n_eff)
    remove(m2_loo)
    
    if (!any(is.na(exp(m2_loglik))))
    {
      m2_rel_n_eff <- relative_eff(exp(m2_loglik))
      suppressWarnings(m2_loo <- loo(m2_loglik, r_eff = m2_rel_n_eff, cores = 4))
    }
    
    c1 <- targetname
    
    if (exists("m1_loo"))
      c2 <- paste0(round(m1_loo$estimates[1,1], 3), " (", round(m1_loo$estimates[1,2], 3), ")")
    else
      c2 <- "NA"
    
    if (exists("m2_loo"))
      c3 <- paste0(round(m2_loo$estimates[1,1], 3), " (", round(m2_loo$estimates[1,2], 3), ")")
    else
      c3 <- "NA"
    
    comparison[n,1] <- c1
    comparison[n,2] <- c2
    comparison[n,3] <- c3
    n <- n + 1
  }
  
  return(comparison)
}

##################################################

mebn.linpred <- function(X, beta, Z, b, g_alpha)
{
  
  mu <- beta %*% X + b %*% Z
  
  # ?? 
  return (dgamma(g_alpha, g_alpha(mu)))
}

##################################################

mebn.AR_comparison <- function(target_variables, graphdir)
{
  library(rstan)
  
  ar_table<-data.frame(matrix(nrow=nrow(target_variables), ncol=4))
  colnames(ar_table) <- c("distribution", "AR1", "CI95l", "CI95u")
  
  n <- 1
  for (targetname in target_variables$Name)
  {
    m1 <- mebn.get_localfit(paste0(graphdir, "/", targetname))
    m1_extract <- extract(m1, pars = c("ar1"))
    
    c2 <- "NA"
    c3 <- "NA"
    c4 <- "NA"

    #if (exists("m1_extract$ar1"))
    #{
      c2 <- mean(m1_extract$ar1)
      # TODO: CI      
    #}
    
    ar_table[n,1] <- targetname
    ar_table[n,2] <- c2
    ar_table[n,3] <- c3
    ar_table[n,4] <- c3
    n <- n + 1
  }
  
  return(ar_table)
}

##################################################

mebn.target_dens_overlays <- function(localfit_directory, target_variables, dataset)
{
  library(rstan)
  library(bayesplot)
  library(ggplot2)
  
  # TODO: Check if dir exists (localfit_directory)
  
  theme_set(bayesplot::theme_default())
  bayesplot::color_scheme_set("purple")
  
  dens_plots <- list()
  i <- 1
  
  for (targetname in target_variables$Name)
  {
    target_blmm <- mebn.get_localfit(paste0(localfit_directory,targetname))
    true_value <- as.vector(dataset[,targetname])
    
    posterior <- extract(target_blmm, pars = c("Y_rep"))
    posterior_y_50 <- posterior$Y_rep[1:50,]
    
    scalemin <- as.numeric(as.character(target_variables[target_variables$Name == targetname,]$ScaleMin))
    scalemax <- as.numeric(as.character(target_variables[target_variables$Name == targetname,]$ScaleMax))
    
    dens_plots[[i]] <- ppc_dens_overlay(true_value, posterior_y_50) + 
      coord_cartesian(xlim = c(scalemin,scalemax)) +
      ggtitle(targetname)
    
    i <- i + 1
  }
  
  bayesplot_grid(plots = dens_plots, legends = FALSE)
}

##################################################

mebn.sampling <- function(inputdata, predictor_columns, target_column, group_column, local_model_cache = "models", stan_mode_file = "BLMM.stan", normalize_values = TRUE, reg_params = NULL)
{
  require(rstan)
  
  # Run Stan parallel on multiple cores
  rstan_options (auto_write=TRUE)
  options (mc.cores=parallel::detectCores ()) 
  
  target_name <- as.vector(target_column$Name)

  localfit <- mebn.get_localfit(target_name, local_model_cache)
  
  # Use cached model if it exists
  if (is.null(localfit))
  {
    stanDat <- mebn.set_model_parameters(predictor_columns, target_column, group_column, inputdata, normalize_values, reg_params)

    initlist <- list(list(beta_Intercept=50),list(beta_Intercept=2),list(beta_Intercept=2),list(beta_Intercept=2))    
    localfit <- stan(file=stan_mode_file, data=stanDat, warmup = 1000, iter=2000, chains=4, init=0, control = list(adapt_delta = 0.80, max_treedepth = 12))
    
    modelcache <- paste0(local_model_cache, "/", target_name, "_blmm", ".rds")
    
    saveRDS(localfit, file=modelcache)
  }
  
  return(localfit)
}  

##################################################

mebn.variational_inference <- function(inputdata, predictor_columns, target_column, group_column, local_model_cache = "models", stan_mode_file = "BLMM.stan", normalize_values = TRUE, reg_params = NULL)
{
  require(rstan)
  
  # Run Stan parallel on multiple cores
  rstan_options (auto_write=TRUE)
  options (mc.cores=parallel::detectCores ()) 
  
  target_name <- as.vector(target_column$Name)
  
  localfit <- mebn.get_localfit(target_name, local_model_cache)
  
  # Use cached model if it exists
  if (is.null(localfit))
  {
    stanDat <- mebn.set_model_parameters(predictor_columns, target_column, group_column, inputdata, normalize_values, reg_params)
    localmodel <- stan_model(file = stan_mode_file)
    localfit <- vb(localmodel, data=stanDat, output_samples=2500, iter=10000)
    
    modelcache <- paste0(local_model_cache, "/", target_name, "_blmm", ".rds")
    saveRDS(localfit, file=modelcache)
  }
  
  return(localfit)
}  

##################################################

mebn.get_rootnodes <- function(g)
{
  which(sapply(sapply(V(g), function(x) neighbors(g,x, mode="in")), length) == 0)
}

##################################################

mebn.predict_graph <- function(reaction_graph, personal_data, alldata, predictor_columns, assumed_targets, local_model_cache = "models", edge_significance_test)
{
  predictor_names <- as.vector(predictor_columns$Name)
  new_input <- subset(personal_data, select = predictor_names)
  
  for (c in 1:dim(assumed_targets)[1])
  {
    target_column <- assumed_targets[c,]
    target_name <- as.vector(target_column$Name)
    
    localfit <- mebn.get_localfit(target_name, local_model_cache)
    
    if (is.null(localfit))
    {
      stop(paste0("Local model for ", target_name, " not found in '", local_model_cache, "'"))
    }
      
    localsummary <- mebn.localsummary(localfit)
    new_response <- subset(personal_data, select = target_name)
  
    # Predict random-effects with BLUP
    ranefs <- mebn.predict_ranefs(localsummary, new_input, new_response)
  
    # Get variances of BLUP, i.e. prediction errors
    localparams <- mebn.set_model_parameters(assumedpredictors, assumedtargets, "notinuse", alldata, TRUE, TRUE)
    err <- mebn.ranef_BLUP_vars(localparams, localsummary)
    ranef_BLUP_vars <- err$var_BLUP
    
    # - Loop through betas for current target
    for (p in 1:length(predictor_names))
    {
      predictor_name <- predictor_names[p]
      
      predicted_coef <- localsummary$fixef[p] + ranefs[p+1]
      prediction_error <- ranef_BLUP_vars[p]
      
      if (edge_significance_test(predicted_coef) == TRUE)
      {
        # Muutettu confband-kaari n?ytt?m??n betan vahvuutta ja CI esitt?? henk.koht. vaihtelua
        
        reaction_graph <- reaction_graph + edge(c(predictor_name, target_name), 
                                                weight = predicted_coef, 
                                                shape   = "confband", 
                                                mean   = predicted_coef,
                                                l95CI  = prediction_error/2,
                                                u95CI  = prediction_error/2)
        
        # Add random-effect for significant predictors
        #reaction_graph <- reaction_graph + vertex(paste0("b_", predictor_name, "_", target_name), label=paste0("b_", predictor_name), color="#AAAAAA", size = 0.5, shape = "disc")
        #reaction_graph <- reaction_graph + vertex(paste0("b_sigma_", predictor_name, "_", target_name), label="b_sigma", color="#AAAAAA", size = localsummary$ranef_sd[p], shape = "disc")
        
        #reaction_graph <- reaction_graph + edge(paste0("b_sigma_", predictor_name, "_", target_name), paste0("b_", predictor_name, "_", target_name), shape = "arrow", weight = localsummary$ranef_sd[p]) 
        #reaction_graph <- reaction_graph + edge(paste0("b_", predictor_name, "_", target_name), target_name, shape = "arrow", weight=1)
      }
      else
      {
        #print(paste0("pruning out edge ", predictor_name, "->", target_name))
      }
    }
    
  } # loop targets
  
  return(reaction_graph)
}

##################################################

mebn.write_gexf <- function(reaction_graph, gexf_path = "C:\\projects\\Responses\\reactions.gexf")
{
  require(rgexf)
  
  # TODO: https://github.com/jacomyal/sigma.js/wiki/Settings
  
  # - hover descriptions for nodes (current vs bounds) 
  # - beta-weights to edges and labels?
  
  MakeRGBA <- function(RGBstring, alpha)
  {
    strtodec <- function(rgb, b, e) { strtoi(paste0("0x", substr(rgb, b, e))) }
    
    RGBA <- data.frame(strtodec(RGBstring, 2, 3), strtodec(RGBstring, 4, 5), strtodec(RGBstring, 6, 7), alpha)  
    colnames(RGBA) <- c("r", "g", "b", "alpha")  
    
    return(RGBA)
  }
  
  graphdata <- get.data.frame(reaction_graph)
  
  nodeviz <- list(color = MakeRGBA(V(reaction_graph)$color, 1.0), size = V(reaction_graph)$size, shape = V(reaction_graph)$shape)
  edgeviz <- list(shape = E(reaction_graph)$shape)
  
  edgesatt <- data.frame(E(reaction_graph)$mean, E(reaction_graph)$l95CI, E(reaction_graph)$u95CI, E(reaction_graph)$b_sigma)
  
  if (length(edgesatt) == 4)
  {
    colnames(edgesatt) <- c("mean", "l95CI", "u95CI", "b_sigma")
  }
  
  nodesatt <- data.frame(V(reaction_graph)$mean, V(reaction_graph)$l95CI, V(reaction_graph)$u95CI, V(reaction_graph)$scalemin, V(reaction_graph)$scalemax)
  if (length(nodesatt) == 5)  
  {
    colnames(nodesatt) <- c("mean", "l95CI", "u95CI", "scalemin", "scalemax")
  }
  
  edgelabels <- data.frame(paste0("beta = ", round(E(reaction_graph)$mean, 3), ", b_sigma = ", round(2*(E(reaction_graph)$u95CI - E(reaction_graph)$mean), 3)))
  
  write.gexf(
    defaultedgetype = "directed",
    nodes = data.frame(V(reaction_graph)$name, V(reaction_graph)$label),
    edges = get.edgelist(reaction_graph),
    edgesWeight = graphdata[,3],
    edgesLabel = edgelabels,
    nodesVizAtt = nodeviz,
    edgesVizAtt = edgeviz,
    edgesAtt = edgesatt,
    nodesAtt = nodesatt,
    output = gexf_path
  )
}

###################################

mebn.BetaLevelTest <- function(LocalModelSummary, PredictorId)
{
  abs(LocalModelSummary$fixef[PredictorId]) > 0.001
}

###################################

mebn.RanefTest <- function(localsummary, PredictorId)
{
  abs(localsummary$fixef[PredictorId]) > 0.001 ||
    abs(localsummary$ranef_sd[PredictorId]) > 0.001
}

###################################

mebn.PersonalSignificanceTest <- function(personal_coef)
{
  abs(personal_coef) > 0.001
}

##################################################

mebn.bipartite_model <- function(reaction_graph, inputdata, predictor_columns, assumed_targets, group_column, local_model_cache, stan_model_file, local_estimation, normalize_values = TRUE, reg_params = NULL)
{
  for (c in 1:dim(assumed_targets)[1])
  {
    target_column <- assumed_targets[c,]
    target_name <- as.vector(target_column$Name)
    
    localfit <- local_estimation(inputdata, predictor_columns, target_column, group_column, local_model_cache, stan_model_file, normalize_values, reg_params)
    
    # Extract model summary
    localsummary <- mebn.localsummary(localfit)

    # - Loop through betas for current target
    predictor_names <- as.vector(predictor_columns$Name)
    
    for (p in 1:length(predictor_names))
    {
      predictor_name <- predictor_names[p]
      
      # Attach the random variable
      reaction_graph <- reaction_graph + edge(c(predictor_name, target_name), 
                                              weight = localsummary$fixef[p], 
                                              b_sigma = localsummary$ranef_sd[p+1],
                                              shape   = "confband", 
                                              mean   = localsummary$fixef[p])
      
      # Fixed-effect
      reaction_graph <- reaction_graph + vertex(paste0("beta_", predictor_name, "_", target_name), 
                                                label=paste0("beta_", predictor_name), 
                                                type="beta", color="#AAAAAA", 
                                                value = localsummary$fixef[p], 
                                                value_l95CI = localsummary$fixef_l95CI[p],
                                                value_u95CI = localsummary$fixef_u95CI[p],
                                                shape = "circle")
      
      reaction_graph <- reaction_graph + edge(paste0("beta_", predictor_name, "_", target_name), paste0("beta_", predictor_name, "_", target_name), shape = "arrow", weight = 1, type = "beta") 
      
      # Add random-effect for significant predictors
      reaction_graph <- reaction_graph + vertex(paste0("b_", predictor_name, "_", target_name), 
                                                label=paste0("b_", predictor_name), 
                                                type="b", 
                                                color="#AAAAAA", 
                                                size = 0.5, 
                                                shape = "circle")
      
      reaction_graph <- reaction_graph + vertex(paste0("b_sigma_", predictor_name, "_", target_name), 
                                                label="b_sigma", 
                                                type="b_sigma", 
                                                color="#AAAAAA", 
                                                value = localsummary$ranef_sd[p],
                                                size = localsummary$ranef_sd[p], 
                                                shape = "circle")
      
      reaction_graph <- reaction_graph + edge(paste0("b_sigma_", predictor_name, "_", target_name), paste0("b_", predictor_name, "_", target_name), shape = "arrow", weight = localsummary$ranef_sd[p], type = "b_sigma") 
      reaction_graph <- reaction_graph + edge(paste0("b_", predictor_name, "_", target_name), target_name, shape = "arrow", weight=1, type = "b")
    }
    
  } # loop targets
  
  return(reaction_graph)
}

##################################################

mebn.personal_graph <- function(person_id, reaction_graph, predictor_columns, assumed_targets, local_model_cache)
{
  library(igraph)
  library(rstan)
  
  for (c in 1:dim(assumed_targets)[1])
  {
    target_column <- assumed_targets[c,]
    target_name <- as.vector(target_column$Name)
    
    localfit <- mebn.get_localfit(paste0(local_model_cache, target_name))

    # extract estimations for specific person
    posterior <- extract(localfit, pars = c("personal_effect", "b"))
    b_blmm <- colMeans(posterior$b)
    beta_b_blmm <- colMeans(posterior$personal_effect)
    
    #pe <- summary(target_blmm, pars="personal_effect", probs = c(0.05, 0.95))$summary[,c(1,4,5)]
    
    # - Loop through betas for current target
    predictor_names <- as.vector(predictor_columns$Name)
    
    for (p in 1:length(predictor_names))
    {
      predictor_name <- predictor_names[p]
      
      # Attach the random variable
      reaction_graph <- reaction_graph + edge(c(predictor_name, target_name),
                                              weight = beta_b_blmm[person_id,p], 
                                              b = b_blmm[person_id,p],
                                              mean = beta_b_blmm[person_id,p])
      
      # Personal effect (beta + b)
      reaction_graph <- reaction_graph + vertex(paste0("personal_", predictor_name, "_", target_name), 
                                                label=paste0("personal_", predictor_name), 
                                                type="personal", color="#AAAAAA", 
                                                value = beta_b_blmm[person_id,p], 
                                                value_l95CI = 0,
                                                value_u95CI = 0,
                                                shape = "circle")
      
      reaction_graph <- reaction_graph + edge(paste0("personal_", predictor_name, "_", target_name), target_name, shape = "arrow", weight = 1, type = "personal") 

      # Personal variation (b)
      reaction_graph <- reaction_graph + vertex(paste0("b_", predictor_name, "_", target_name), 
                                                label=paste0("b_", predictor_name), 
                                                type="b", color="#AAAAAA", 
                                                value = b_blmm[person_id,p], 
                                                value_l95CI = 0,
                                                value_u95CI = 0,
                                                shape = "circle")
      
      reaction_graph <- reaction_graph + edge(paste0("b_", predictor_name, "_", target_name), target_name, shape = "arrow", weight = 1, type = "b") 
    }
    
  } # loop targets
  
  return(reaction_graph)
}

##################################################

mebn.plot_typical_effects <- function(reaction_graph, top_effects)
{
  # Parameter and hyperparameter nodes are removed and visualized otherwise
  visual_graph <- reaction_graph
  
  # Remove edges to latent variables
  visual_graph <- delete.edges(visual_graph, which(E(visual_graph)$type=="beta"))
  visual_graph <- delete.edges(visual_graph, which(E(visual_graph)$type=="b_sigma"))
  visual_graph <- delete.edges(visual_graph, which(E(visual_graph)$type=="b"))
  
  # Remove nodes of latent variable
  visual_graph <- delete.vertices(visual_graph, which(V(visual_graph)$type=="beta"))
  visual_graph <- delete.vertices(visual_graph, which(V(visual_graph)$type=="b_sigma"))
  visual_graph <- delete.vertices(visual_graph, which(V(visual_graph)$type=="b"))
  
  # Filter only the most significant edges having large typical effect or large personal variance
  alledges <- E(visual_graph)
  top_neg_edges <- head(alledges[order(alledges$weight)], top_effects)
  top_pos_edges <- head(alledges[order(-alledges$weight)], top_effects)
  top_pers_edges <- head(alledges[order(-alledges$b_sigma)], top_effects)
  
  # Comment out this row to see all the connections at the model
  visual_graph <- delete.edges(visual_graph, alledges[-c(top_neg_edges, top_pos_edges, top_pers_edges)])
  
  # Graph layout
  V(visual_graph)$size = 10 
  # - put all blood test values in own rank
  bipa_layout <- layout_as_bipartite(visual_graph, types = V(visual_graph)$type == "100")
  # - flip layout sideways, from left to right
  gap <- 4
  bipa_layout <- cbind(bipa_layout[,2]*gap, bipa_layout[,1])

  # Color and size encoding for edges according to beta coefficient
  E(visual_graph)[E(visual_graph)$weight > 0]$color="red"
  E(visual_graph)[E(visual_graph)$weight < 0]$color="blue"
  E(visual_graph)$width = abs(E(visual_graph)$weight) * 7
  
  # Idea: make own edge for personal variance to visualization graph
  
  plot(visual_graph, 
       layout=bipa_layout, 
       rescale=TRUE,
       vertex.label.family="Helvetica",
       vertex.label.color="black",
       vertex.label.cex=1,
       edge.arrow.size=1,
       edge.arrow.width=1)
}

##################################################

mebn.plot_personal_effects <- function(personal_graph, top_effects)
{
  # Parameter and hyperparameter nodes are removed and visualized otherwise
  visual_graph <- personal_graph
  
  # Remove edges to latent variables
  visual_graph <- delete.edges(visual_graph, which(E(visual_graph)$type=="beta"))
  visual_graph <- delete.edges(visual_graph, which(E(visual_graph)$type=="b_sigma"))
  visual_graph <- delete.edges(visual_graph, which(E(visual_graph)$type=="b"))
  visual_graph <- delete.edges(visual_graph, which(E(visual_graph)$type=="personal"))
  
  # Remove nodes of latent variable
  visual_graph <- delete.vertices(visual_graph, which(V(visual_graph)$type=="beta"))
  visual_graph <- delete.vertices(visual_graph, which(V(visual_graph)$type=="b_sigma"))
  visual_graph <- delete.vertices(visual_graph, which(V(visual_graph)$type=="b"))
  visual_graph <- delete.vertices(visual_graph, which(V(visual_graph)$type=="personal"))
  
  # Filter only the most significant edges having large typical effect
  alledges <- E(visual_graph)
  top_neg_edges <- head(alledges[order(alledges$weight)], top_effects)
  top_pos_edges <- head(alledges[order(-alledges$weight)], top_effects)

  # Comment out this row to see all the connections at the model
  visual_graph <- delete.edges(visual_graph, alledges[-c(top_neg_edges, top_pos_edges)])
  
  # Graph layout
  V(visual_graph)$size = 10 
  # - put all blood test values in own rank
  bipa_layout <- layout_as_bipartite(visual_graph, types = V(visual_graph)$type == "100")
  # - flip layout sideways, from left to right
  gap <- 4
  bipa_layout <- cbind(bipa_layout[,2]*gap, bipa_layout[,1])
  
  # Color and size encoding for edges according to beta coefficient
  E(visual_graph)[E(visual_graph)$weight > 0]$color="red"
  E(visual_graph)[E(visual_graph)$weight < 0]$color="blue"
  E(visual_graph)$width = abs(E(visual_graph)$weight) * 7
  
  plot(visual_graph, 
       layout=bipa_layout, 
       rescale=TRUE,
       vertex.label.family="Helvetica",
       vertex.label.color="black",
       vertex.label.cex=1,
       edge.arrow.size=1,
       edge.arrow.width=1)
}

##################################################

mebn.plot_personal_variations <- function(reaction_graph, top_effects)
{
  # Parameter and hyperparameter nodes are removed and visualized otherwise
  visual_graph <- reaction_graph
  
  # Remove edges to latent variables
  visual_graph <- delete.edges(visual_graph, which(E(visual_graph)$type=="beta"))
  visual_graph <- delete.edges(visual_graph, which(E(visual_graph)$type=="b_sigma"))
  visual_graph <- delete.edges(visual_graph, which(E(visual_graph)$type=="b"))
  
  # Remove nodes of latent variable
  visual_graph <- delete.vertices(visual_graph, which(V(visual_graph)$type=="beta"))
  visual_graph <- delete.vertices(visual_graph, which(V(visual_graph)$type=="b_sigma"))
  visual_graph <- delete.vertices(visual_graph, which(V(visual_graph)$type=="b"))
  
  # Filter only the most significant edges having large typical effect or large personal variance
  alledges <- E(visual_graph)
  top_neg_edges <- head(alledges[order(alledges$weight)], top_effects)
  top_pos_edges <- head(alledges[order(-alledges$weight)], top_effects)
  top_pers_edges <- head(alledges[order(-alledges$b_sigma)], top_effects)
  
  # Comment out this row to see all the connections at the model
  visual_graph <- delete.edges(visual_graph, alledges[-c(top_neg_edges, top_pos_edges, top_pers_edges)])
  
  # Graph layout
  V(visual_graph)$size = 10 
  # - put all blood test values in own rank
  bipa_layout <- layout_as_bipartite(visual_graph, types = V(visual_graph)$type == "100")
  # - flip layout sideways, from left to right
  gap <- 4
  bipa_layout <- cbind(bipa_layout[,2]*gap, bipa_layout[,1])
  
  # Color and size encoding for edges according to beta coefficient
  E(visual_graph)$color="gray"
  E(visual_graph)$width = abs(E(visual_graph)$b_sigma) * 7
  
  plot(visual_graph, 
       layout=bipa_layout, 
       rescale=TRUE,
       vertex.label.family="Helvetica",
       vertex.label.color="black",
       vertex.label.cex=1,
       edge.arrow.size=1,
       edge.arrow.width=1)
}

##################################################

mebn.set_evidence <- function(reaction_graph, evidence)
{
  # Values of nodes (like optimal blood test values) can be set to network
  
}

##################################################

mebn.evaluate <- function(reaction_graph, targets)
{
  # After the evidence values have been set, the rest of the network can be evaluated 
  
  # Stan sampling?
  
}

##################################################
