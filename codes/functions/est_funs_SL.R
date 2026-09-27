library(SuperLearner)


SL_est_density_ratio_1vs0 <- function(SL, variable, data1, data0, sd = 1){
  
  density_1 = dnorm(variable, mean = predict(SL, data1, onlySL = T)$pred, sd = sd)
  density_0 = dnorm(variable, mean = predict(SL, data0, onlySL = T)$pred, sd = sd)
  
  return(density_1/density_0)
}



SL_est_bayes_ratio_1vs0 <- function(SL, data1, data0){
  
  prob1 <- predict(SL, data1, onlySL = T)$pred       # T1 = 1; or (T1 = 1 & T2 = 1)
  prob0 <- 1 - predict(SL, data0, onlySL = T)$pred   # T1 = 0; or (T1 = 0 & T2 = 0)
  
  return(prob1/prob0)
}


SL_est_ratio_1vs0 <- function(SL, predict_data){
  
  prob1 <- predict(SL, predict_data, onlySL = T)$pred       # T = 1| ....
  prob0 <- 1 - prob1                                        # T = 0| ....
  
  return(prob1/prob0)
}


SL_est_outcome_regressions <- function(Y = "Y", data, pred.data.list, predictor.list, 
                                       SL_library = SL_library) {
  
  # Initialize list to store predictions
  outcome_pred <- vector("list", length(predictor.list) + 1)
  
  # Store original Y values
  outcome_pred[[1]] <- pred_Y <- data[[Y]]
  
  # Loop through predictor sets and fit models using SuperLearner
  for (i in seq_along(predictor.list)) {
    
    predictors <- predictor.list[[i]]
  
    model <- SuperLearner(Y = pred_Y, X = data[, predictors, drop = FALSE], 
                          family = gaussian(), SL.library = SL_library)
  
    pred_Y <- predict(model, newdata = pred.data.list[[i]][, predictors, drop = FALSE], onlySL = T)$pred
    outcome_pred[[i + 1]] <- pred_Y
  }
  
  # Assign names to outcome_pred list
  names(outcome_pred) <- c('Y', 'u', paste0("B", seq_along(predictor.list)[-1]-1))
  
  return(outcome_pred)
}




# *******************************************************
#  SuperLearner estimation        ----------------------- 
# *******************************************************

SL_nuisance_funs <-function(data, X = c(),
                             density = T, bayes = T, joint = "OR", 
                             SL_library = c('SL.glm', 'SL.glm.interaction')){
  
  # *******************************************************
  # 0. prepare data                 ----------------------- 
  
  head(data)
  
  data1 <- data |> mutate(T1 = 1, T2 = 1)
  data0 <- data |> mutate(T1 = 0, T2 = 0)
  
  cum_var <- Reduce(c, c(list(X), 'T1', 'M1', 'L1', 'T2', 'M2', 'L2'), accumulate = TRUE)
  names(cum_var) = c('c.X', 'c.T1', 'c.M1', 'c.L1', 'c.T2', 'c.M2', 'c.L2')
  
  cum_var_noT1 <- Reduce(c, c(list(X), 'M1', 'L1', 'T2', 'M2', 'L2'), accumulate = TRUE)
  names(cum_var_noT1) = c('c.X',  'c.M1', 'c.L1', 'c.T2', 'c.M2', 'c.L2')
  
  cum_var_noT2 <- Reduce(c, c(list(X), 'T1', 'M1', 'L1', 'M2', 'L2'), accumulate = TRUE)
  names(cum_var_noT2) = c('c.X', 'c.T1', 'c.M1', 'c.L1', 'c.M2', 'c.L2')
  
  cum_var_noT <- Reduce(c, c(list(X), 'M1', 'L1', 'M2', 'L2'), accumulate = TRUE)
  names(cum_var_noT) = c('c.X', 'c.M1', 'c.L1', 'c.M2', 'c.L2')
  
  
  
  # *******************************************************
  # 1. nuisance functions           ----------------------- 
  
  ## (1) propensity functions ------------
  pro.T1 = SuperLearner(Y = data$T1, X = data[ , X], 
                        family = binomial(), SL.library = SL_library)
  pro.T2 = SuperLearner(Y = data$T2, X = data[ , cum_var$c.L1], 
                        family = binomial(), SL.library = SL_library)
  
  prob1.T1 = pro.T1$SL.predict
  prob0.T1 = 1 - prob1.T1
  
  prob1.T2 = predict(pro.T2, data1[ , cum_var$c.L1], onlySL = T)$pred        # T2 = 1 & T1 = 1
  prob0.T2 = 1 - predict(pro.T2, data0[ , cum_var$c.L1], onlySL = T)$pred    # T2 = 0 & T1 = 0
  
  
  ## (2) density ------------
  if(density){
    
    f_M1 = SuperLearner(Y = data$M1, X = data[ , cum_var$c.T1], 
                        family = gaussian(), SL.library = SL_library)
    f_L1 = SuperLearner(Y = data$L1, X = data[ , cum_var$c.M1], 
                        family = gaussian(), SL.library = SL_library)
    f_M2 = SuperLearner(Y = data$M2, X = data[ , cum_var$c.T2], 
                        family = gaussian(), SL.library = SL_library)
    f_L2 = SuperLearner(Y = data$L2, X = data[ , cum_var$c.M2], 
                        family = gaussian(), SL.library = SL_library)
    
    
    density.ratio.M1 = SL_est_density_ratio_1vs0(f_M1, data$M1, data1[,cum_var$c.T1], data0[,cum_var$c.T1])
    density.ratio.L1 = SL_est_density_ratio_1vs0(f_L1, data$L1, data1[,cum_var$c.M1], data0[,cum_var$c.M1])
    density.ratio.M2 = SL_est_density_ratio_1vs0(f_M2, data$M2, data1[,cum_var$c.T2], data0[,cum_var$c.T2])
    density.ratio.L2 = SL_est_density_ratio_1vs0(f_L2, data$L2, data1[,cum_var$c.M2], data0[,cum_var$c.M2])
    
  }else{
    density.ratio.M1 = density.ratio.L1 = density.ratio.M2 = density.ratio.L2 = NULL
  }
  
  ## (3) bayes ------------

  if(bayes == FALSE){
    
    bayes.ratio.T1_M1 = bayes.ratio.T1_L1 = bayes.ratio.T1_M2 = bayes.ratio.T1_L2  = 
      bayes.ratio.T2_M2 = bayes.ratio.T2_L2 = NULL
    
  }else if(joint == "chain"){
    
    f_T1_M1 <- SuperLearner(Y = data$T1, X = data[,cum_var_noT$c.M1], 
                            family = binomial(), SL.library = SL_library)
    f_T1_L1 <- SuperLearner(Y = data$T1, X = data[,cum_var_noT$c.L1], 
                            family = binomial(), SL.library = SL_library)
    f_T1_M2 <- SuperLearner(Y = data$T1, X = data[,cum_var_noT$c.M2], 
                            family = binomial(), SL.library = SL_library)
    f_T1_L2 <- SuperLearner(Y = data$T1, X = data[,cum_var_noT$c.L2], 
                            family = binomial(), SL.library = SL_library)
    
    f_T2_M2 <- SuperLearner(Y = data$T2, X = data[,cum_var_noT2$c.M2], 
                            family = binomial(), SL.library = SL_library)
    f_T2_L2 <- SuperLearner(Y = data$T2, X = data[,cum_var_noT2$c.L2], 
                            family = binomial(), SL.library = SL_library)
    
    bayes.ratio.T1_M1 = SL_est_ratio_1vs0(f_T1_M1, data[,cum_var_noT$c.M1])
    bayes.ratio.T1_L1 = SL_est_ratio_1vs0(f_T1_L1, data[,cum_var_noT$c.L1])
    bayes.ratio.T1_M2 = SL_est_ratio_1vs0(f_T1_M2, data[,cum_var_noT$c.M2])
    bayes.ratio.T1_L2 = SL_est_ratio_1vs0(f_T1_L2, data[,cum_var_noT$c.L2])
    
    # (t1 = t2 = 0 or 1)
    bayes.ratio.T2_M2 = SL_est_bayes_ratio_1vs0(f_T2_M2, data1[,cum_var_noT2$c.M2], data0[,cum_var_noT2$c.M2])
    bayes.ratio.T2_L2 = SL_est_bayes_ratio_1vs0(f_T2_L2, data1[,cum_var_noT2$c.L2], data0[,cum_var_noT2$c.L2])    
    
  }else if(joint == "OR"){
    
    f_T1_M1 <- SuperLearner(Y = data$T1, X = data[,cum_var_noT$c.M1], 
                            family = binomial(), SL.library = SL_library)
    f_T1_L1 <- SuperLearner(Y = data$T1, X = data[,cum_var_noT$c.L1], 
                            family = binomial(), SL.library = SL_library)
    
    f_T1_M2 <- SuperLearner(Y = data$T1, X = data[,cum_var_noT1$c.M2], 
                            family = binomial(), SL.library = SL_library)
    f_T1_L2 <- SuperLearner(Y = data$T1, X = data[,cum_var_noT1$c.L2], 
                            family = binomial(), SL.library = SL_library)
    
    f_T2_M2 <- SuperLearner(Y = data$T2, X = data[,cum_var_noT2$c.M2], 
                            family = binomial(), SL.library = SL_library)
    f_T2_L2 <- SuperLearner(Y = data$T2, X = data[,cum_var_noT2$c.L2], 
                            family = binomial(), SL.library = SL_library)
    
    bayes.ratio.T1_M1 = SL_est_ratio_1vs0(f_T1_M1, predict_data = data[,cum_var_noT$c.M1])
    bayes.ratio.T1_L1 = SL_est_ratio_1vs0(f_T1_L1, predict_data = data[,cum_var_noT$c.L1])
    
    # T1 | M1.L1.T2=0.M2
    bayes.ratio.T1_M2 = SL_est_ratio_1vs0(f_T1_M2, predict_data = data0[,cum_var_noT1$c.M2])
    # T1 | M1.L1.T2=0.M2.L2
    bayes.ratio.T1_L2 = SL_est_ratio_1vs0(f_T1_L2, predict_data = data0[,cum_var_noT1$c.L2])
    
    # T2 | T1=1.M1.L1.M2 
    bayes.ratio.T2_M2 = SL_est_ratio_1vs0(f_T2_M2, predict_data = data1[,cum_var_noT2$c.M2])
    # T2 | T1=1.M1.L1.M2.L2
    bayes.ratio.T2_L2 = SL_est_ratio_1vs0(f_T2_L2, predict_data = data1[,cum_var_noT2$c.L2])       
  }
  
  # (4) outcome regression ------------
  pred.data.list <- list(
    data0, data1, data0, data1, data0
  )
  
  # use predY
  predictor.list <- list(
    cum_var$c.L2,
    cum_var$c.M2,
    cum_var$c.T2,
    cum_var$c.M1,
    cum_var$c.T1
  )
  
 
  out_pred <- SL_est_outcome_regressions(Y = "Y", data, pred.data.list, predictor.list, 
                                         SL_library = SL_library) 
    
  
  
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