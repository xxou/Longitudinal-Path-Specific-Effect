# *****************************************************************
# One step estimator for longitudinal treatments
# *****************************************************************
# DGP:  Feb - 22 - 2025 Simplest data
# Xiaxian Ou
# two treatments T1 -> M1 -> L1 -> T2 -> M2 -> L2 -> Y
# *****************************************************************

library(tidyverse)

lm_est_density_ratio_1vs0 <- function(mu_formula, variable, data, data1, data0){
  fit <- lm(mu_formula, data)
  
  density_1 = dnorm(data[[variable]], mean = predict(fit, data1), sd = summary(fit)$sigma)
  density_0 = dnorm(data[[variable]], mean = predict(fit, data0), sd = summary(fit)$sigma)
  
  return(density_1/density_0)
}



# bayes rule + chain rule 
glm_est_bayes_ratio_1vs0 <- function(T_formula, data, data1, data0){
  fit <- glm(T_formula, data, family = binomial())
  
  prob1 <- predict(fit, data1, type = "response")       # T1 = 1; or (T1 = 1 & T2 = 1)
  prob0 <- 1 - predict(fit, data0, type = "response")   # T1 = 0; or (T1 = 0 & T2 = 0)
  
  return(prob1/prob0)
}



# joint OR
# (prob)/(1-prob)
glm_est_ratio_1vs0 <- function(T_formula, data, predict_data, At_correction = FALSE){
  if(At_correction){
    fit <- glm(T_formula, data, family = binomial(), offset = At)
  }else{
    fit <- glm(T_formula, data, family = binomial())
  }
  
  
  prob1 <- predict(fit, predict_data, type = "response")  
  prob0 <- 1 - prob1   
  
  return(prob1/prob0)
}



  


est_outcome_regressions <- function(Y = "Y", data, pred.data.list, formulas.list) {
  # Initialize list to store predictions
  outcome_pred <- vector("list", length(formulas.list) + 1)
  
  # Store original Y values
  outcome_pred[[1]] <- data[[Y]]
  data$predY <- data[[Y]]
  
  # Loop through formulas and fit models
  for (i in seq_along(formulas.list)) {
    model <- glm(formula = formulas.list[[i]], data, family = gaussian())  # Fit model
    outcome_pred[[i + 1]]  <- predict(model, pred.data.list[[i]], type = "response")  # Predict on corresponding dataset
    data$predY <- outcome_pred[[i + 1]]
  }
  
  names(outcome_pred) = c('Y', names(formulas.list))
  return(outcome_pred)
}



GLM_nuisance_funs <-function(data, f.pro.T.list, f.mediator.list = NULL, f.bayes.list = NULL, seq.reg.list, joint = "OR", At_correction = FALSE){
  
  
  list2env(f.pro.T.list, envir = environment())
  if(!is.null(f.mediator.list)){list2env(f.mediator.list, envir = environment())} 
  if(!is.null(f.bayes.list)){list2env(f.bayes.list, envir = environment())} 
  
  # *******************************************************
  # 0. prepare data                 ----------------------- 
  
  head(data)
  
  data1 <- data |> mutate(T1 = 1, T2 = 1)
  data0 <- data |> mutate(T1 = 0, T2 = 0)
  
  # *******************************************************
  # 1. nuisance functions           ----------------------- 

  ## (1) propensity functions ------------
  pro.T1 = glm(pro.T1_formula, data, family = binomial())
  pro.T2 = glm(pro.T2_formula, data, family = binomial())
  
  prob1.T1 = predict(pro.T1,  type = "response")           # T1 = 1 | X
  prob0.T1 = 1 - prob1.T1                                  # T1 = 0 | X
  
  prob1.T2 = predict(pro.T2, data1, type = "response")      # T2 = 1 | T1 = 1, L1, M1, X
  prob0.T2 = 1 - predict(pro.T2, data0, type = "response")  # T2 = 0 | T1 = 0, L1, M1, X
 
  
  ## (2) density ------------
  if(is.null(f.mediator.list)){
    
    density.ratio.M1 = density.ratio.L1 = density.ratio.M2 = density.ratio.L2 = NULL
    
  }else{
    
    density.ratio.M1 = lm_est_density_ratio_1vs0(M1_formula, "M1", data, data1, data0)
    density.ratio.L1 = lm_est_density_ratio_1vs0(L1_formula, "L1", data, data1, data0)
    density.ratio.M2 = lm_est_density_ratio_1vs0(M2_formula, "M2", data, data1, data0)
    density.ratio.L2 = lm_est_density_ratio_1vs0(L2_formula, "L2", data, data1, data0)
  
  }
  
  
  ## (3) bayes ------------
  if(is.null(f.bayes.list)){
    
    bayes.ratio.T1_M1 = bayes.ratio.T1_L1 = bayes.ratio.T1_M2 = bayes.ratio.T1_L2  = 
      bayes.ratio.T2_M2 = bayes.ratio.T2_L2 = NULL
    
  }else if(joint == "chain"){
    
    
    bayes.ratio.T1_M1 = glm_est_bayes_ratio_1vs0(T1_M1_formula, data, data1, data0)
    bayes.ratio.T1_L1 = glm_est_bayes_ratio_1vs0(T1_L1_formula, data, data1, data0)   # T1 | M1.L1
    bayes.ratio.T1_M2 = glm_est_bayes_ratio_1vs0(T1_M2_formula, data, data1, data0)   # T1 | M1.L1.M2 (no T2)
    bayes.ratio.T1_L2 = glm_est_bayes_ratio_1vs0(T1_L2_formula, data, data1, data0)   # T1 | M1.L1.M2.L2 (no T2)
    
    
    # (t1 = t2 = 0 or 1)
    bayes.ratio.T2_M2 = glm_est_bayes_ratio_1vs0(T2_M2_formula, data, data1, data0)   # T2 | T1.M1.L1.M2 
    bayes.ratio.T2_L2 = glm_est_bayes_ratio_1vs0(T2_L2_formula, data, data1, data0)   # T2 | T1.M1.L1.M2.L2
    
    
  }else if(joint == "OR"){
    
    bayes.ratio.T1_M1 = glm_est_ratio_1vs0(T1_M1_formula, data, predict_data = data, At_correction = FALSE)     # T1 | M1
    bayes.ratio.T1_L1 = glm_est_ratio_1vs0(T1_L1_formula, data, predict_data = data, At_correction = FALSE)     # T1 | M1.L1
    
   
    # At = -log [ p(T2 = t | T1 = 0, L1, M1, X)/p(T2 = t | T1 = 1, L1, M1, X) ]
    # prob0.T2: p(T2 = 0 | T1 = 0, L1, M1, X)
    # prob1.T2: p(T2 = 1 | T1 = 1, L1, M1, X)
    data$At = -log( ( data$T2*(1-prob0.T2) + (1-data$T2)*prob0.T2 )/( data$T2*prob1.T2 + (1-data$T2)*(1-prob1.T2) ) )
    # At=0 = -log [ p(T2 = 0 | T1 = 0, L1, M1, X)/p(T2 = 0 | T1 = 1, L1, M1, X) ]
    data0$At = -log( prob0.T2 / (1-prob1.T2) )
    
    
    bayes.ratio.T1_M2 = glm_est_ratio_1vs0(T1_M2_formula, data, predict_data = data0, At_correction = At_correction)    # T1 | M1.L1.T2=0.M2 
    bayes.ratio.T1_L2 = glm_est_ratio_1vs0(T1_L2_formula, data, predict_data = data0, At_correction = At_correction)    # T1 | M1.L1.T2=0.M2.L2
    
  
    bayes.ratio.T2_M2 = glm_est_ratio_1vs0(T2_M2_formula, data, predict_data = data1, At_correction = FALSE)   # T2 | T1=1.M1.L1.M2 
    bayes.ratio.T2_L2 = glm_est_ratio_1vs0(T2_L2_formula, data, predict_data = data1, At_correction = FALSE)   # T2 | T1=1.M1.L1.M2.L2
  
  }
  
  # (4) outcome regression ------------
  pred.data.list <- list(
    data0, data1, data0, data1, data0
  )
  
  out_pred <- est_outcome_regressions(Y = "Y", data, pred.data.list, formulas.list = seq.reg.list)
  
  
  # *******************************************************
  # 3. output                 ----------------------- 
  
  out <- list(
    prob1.T1 = prob1.T1, prob0.T1 = prob0.T1, 
    prob1.T2 = prob1.T2, prob0.T2 = prob0.T2,
    
    density.ratio.M1 = density.ratio.M1,
    density.ratio.L1 = density.ratio.L1,
    density.ratio.M2 = density.ratio.M2,
    density.ratio.L2 = density.ratio.L2,
    
    bayes.ratio.T1_M1 = bayes.ratio.T1_M1,
    bayes.ratio.T1_L1 = bayes.ratio.T1_L1,
    bayes.ratio.T1_M2 = bayes.ratio.T1_M2,
    bayes.ratio.T1_L2 = bayes.ratio.T1_L2,
    
    bayes.ratio.T2_M2 = bayes.ratio.T2_M2,
    bayes.ratio.T2_L2 = bayes.ratio.T2_L2,
    
    Y = out_pred$Y,
    u = out_pred$u,
    B1 = out_pred$B1,
    B2 = out_pred$B2,
    B3 = out_pred$B3,
    B4 = out_pred$B4
  )
  
  
  return(out)
}
