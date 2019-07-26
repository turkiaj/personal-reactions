
# Helper fuctions for model evaluation
library(dplyr)
library(tidyr)

true_values <- sysdimet[c("SUBJECT_ID", "WEEK", as.vector(assumedtargets$Name))]

get_true_response <- function(holdout_number) {
  
  #holdout_subject <- sprintf("S%02d", holdout_number)
  holdout_subject <- levels(sysdimet$SUBJECT_ID)[holdout_number]
  
  true_response <- true_values %>%
    filter(SUBJECT_ID == holdout_subject) %>%
    select(as.vector(assumedtargets$Name)) %>%
    t %>%
    as.matrix
  
  colnames(true_response) <- c(0,4,8,12)
  return(true_response)
}

get_pred_response <- function(localfit_directory) {
  
  pred_response <- NULL
  
  for (targetname in assumedtargets$Name)
  {
    print(targetname)
    target_blmm <- mebn.get_localfit(paste0(localfit_directory,targetname))
    ms <- rstan::summary(target_blmm, pars=c("Y_pred"), probs=c(0.10, 0.90), na.rm = TRUE)
    
    if (is.null(pred_response))
      pred_response <- as.vector(ms$summary[1:4,c(1)])
    else
      pred_response <- rbind(pred_response, as.vector(ms$summary[1:4,c(1)]))
  }
  
  colnames(pred_response) <- c(0,4,8,12)
  rownames(pred_response) <- as.vector(assumedtargets$Name)
  
  return(pred_response)
}

get_rep_response_for_all <- function(localfit_directory) {
  
  rep_response <- NULL
  
  for (targetname in assumedtargets$Name)
  {
    target_blmm <- mebn.get_localfit(paste0(localfit_directory,targetname))
    ms <- rstan::summary(target_blmm, pars=c("Y_rep"), probs=c(0.10, 0.90), na.rm = TRUE)
    
    if (is.null(rep_response))
      rep_response <- as.vector(ms$summary[,c(1)])
    else
      rep_response <- rbind(rep_response, as.vector(ms$summary[,c(1)]))
  }
  
  rownames(rep_response) <- as.vector(assumedtargets$Name)
  
  return(rep_response)
}

get_rep_response <- function(rep_matrix, subject_number) {
  
  s <- (subject_number-1)*4 + 1
  e <- s + 3
  
  return(rep_matrix[,s:e])
}  

make_delta <- function(true_response)
{
  true_delta <- true_response
  true_delta[,4] <- true_delta[,4] - true_delta[,3]
  true_delta[,3] <- true_delta[,3] - true_delta[,2]
  true_delta[,2] <- true_delta[,2] - true_delta[,1]
  true_delta[,1] <- 0
  
  return(true_delta)
}