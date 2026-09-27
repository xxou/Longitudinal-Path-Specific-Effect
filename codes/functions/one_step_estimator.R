# *****************************************************************
# One step estimator for longitudinal treatments
# *****************************************************************
# Xiaxian Ou
# two treatments T1 -> M1 -> L1 -> T2 -> M2 -> L2 -> Y
# *****************************************************************


condtion.funcs <- function(true_funcs, false_funcs, type = "bayes"){
  
  
  prob_names <- c("prob1.T1", "prob0.T1", "prob1.T2", "prob0.T2")
  
  false_funcs <- purrr::map_if(
    false_funcs, names(false_funcs) %in% prob_names,  ~ pmax(pmin(.x, 0.99), 0.01)  
    )
  
 
  # *****************************
  ##   density   ----------------
  # *****************************
  
  # density ratio
  # condition 1. f(T1|..), f(T2|..), f(L1|..), f(L2|..)
  # condition 2. outcome regression E(Y|..), f(M1|..), f(M2|..), f(L1|..), f(L2|..)
  # condition 3. f(T1|..), f(T2|..), E(Y|..), f(M1|..), f(M2|..)
  
  if(type == "density"){
    d.con1_funs <- d.con2_funs <- d.con3_funs <- false_funcs
    
    # ******* condition 1 ***********
    d.con1_funs$prob1.T1 <- true_funcs$prob1.T1; d.con1_funs$prob0.T1 <- true_funcs$prob0.T1
    d.con1_funs$prob1.T2 <- true_funcs$prob1.T2; d.con1_funs$prob0.T2 <- true_funcs$prob0.T2
    
    d.con1_funs$density.ratio.L1 <- true_funcs$density.ratio.L1
    d.con1_funs$density.ratio.L2 <- true_funcs$density.ratio.L2
    
    
    # ******* condition 2 ***********
    d.con2_funs$u <- true_funcs$u
    d.con2_funs$B1 <- true_funcs$B1
    d.con2_funs$B2 <- true_funcs$B2
    d.con2_funs$B3 <- true_funcs$B3
    d.con2_funs$B4 <- true_funcs$B4
    
    d.con2_funs$density.ratio.M1 <- true_funcs$density.ratio.M1
    d.con2_funs$density.ratio.M2 <- true_funcs$density.ratio.M2
    
    d.con2_funs$density.ratio.L1 <- true_funcs$density.ratio.L1
    d.con2_funs$density.ratio.L2 <- true_funcs$density.ratio.L2
    
    
    # ******* condition 3 ***********
    d.con3_funs$prob1.T1 <- true_funcs$prob1.T1; d.con3_funs$prob0.T1 <- true_funcs$prob0.T1
    d.con3_funs$prob1.T2 <- true_funcs$prob1.T2; d.con3_funs$prob0.T2 <- true_funcs$prob0.T2
    
    d.con3_funs$u <- true_funcs$u
    d.con3_funs$density.ratio.M1 <- true_funcs$density.ratio.M1
    d.con3_funs$density.ratio.M2 <- true_funcs$density.ratio.M2
    
    # ******* output ***********
    out <- list(d.con1_funs = d.con1_funs, d.con2_funs = d.con2_funs, d.con3_funs = d.con3_funs)
  }
  
  
  # *****************************
  ##     bayes   ----------------
  # *****************************
  
  if(type == "bayes"){
    
    b.con1_funs <- b.con2_funs <- b.con3_funs <- b.con4_funs <- b.con5_funs <- b.con6_funs <- false_funcs
    
    
    # ******* condition 1 ***********
    # EIF line 1
    b.con1_funs$prob1.T1 <- true_funcs$prob1.T1; b.con1_funs$prob0.T1 <- true_funcs$prob0.T1
    b.con1_funs$prob1.T2 <- true_funcs$prob1.T2; b.con1_funs$prob0.T2 <- true_funcs$prob0.T2
    
    b.con1_funs$bayes.ratio.T1_M1 = true_funcs$bayes.ratio.T1_M1
    b.con1_funs$bayes.ratio.T1_L1 = true_funcs$bayes.ratio.T1_L1
    b.con1_funs$bayes.ratio.T1_M2 = true_funcs$bayes.ratio.T1_M2
    b.con1_funs$bayes.ratio.T1_L2 = true_funcs$bayes.ratio.T1_L2
    
    b.con1_funs$bayes.ratio.T2_L1 = true_funcs$bayes.ratio.T2_L1
    b.con1_funs$bayes.ratio.T2_M2 = true_funcs$bayes.ratio.T2_M2
    b.con1_funs$bayes.ratio.T2_L2 = true_funcs$bayes.ratio.T2_L2
    
    
    # ******* condition 2 ***********
    # EIF line 2
    
    b.con2_funs$u <- true_funcs$u
    
    b.con2_funs$prob1.T1 <- true_funcs$prob1.T1; b.con2_funs$prob0.T1 <- true_funcs$prob0.T1
    b.con2_funs$prob1.T2 <- true_funcs$prob1.T2; b.con2_funs$prob0.T2 <- true_funcs$prob0.T2
    
    b.con2_funs$bayes.ratio.T1_M1 = true_funcs$bayes.ratio.T1_M1
    b.con2_funs$bayes.ratio.T1_L1 = true_funcs$bayes.ratio.T1_L1
    b.con2_funs$bayes.ratio.T1_M2 = true_funcs$bayes.ratio.T1_M2
    
    b.con2_funs$bayes.ratio.T2_L1 = true_funcs$bayes.ratio.T2_L1
    b.con2_funs$bayes.ratio.T2_M2 = true_funcs$bayes.ratio.T2_M2
    
    
    
    # ******* condition 3 ***********
    # EIF line 3
    
    b.con3_funs$u <- true_funcs$u
    b.con3_funs$B1 <- true_funcs$B1
    
    b.con3_funs$prob1.T1 <- true_funcs$prob1.T1; b.con3_funs$prob0.T1 <- true_funcs$prob0.T1
    b.con3_funs$prob1.T2 <- true_funcs$prob1.T2; b.con3_funs$prob0.T2 <- true_funcs$prob0.T2
    
    b.con3_funs$bayes.ratio.T1_M1 = true_funcs$bayes.ratio.T1_M1
    b.con3_funs$bayes.ratio.T1_L1 = true_funcs$bayes.ratio.T1_L1
    

    # ******* condition 4 ***********
    # EIF line 4
    
    b.con4_funs$u <- true_funcs$u
    b.con4_funs$B1 <- true_funcs$B1
    b.con4_funs$B2 <- true_funcs$B2
    
    b.con4_funs$prob1.T1 <- true_funcs$prob1.T1; b.con4_funs$prob0.T1 <- true_funcs$prob0.T1
    
    b.con4_funs$bayes.ratio.T1_M1 = true_funcs$bayes.ratio.T1_M1
    
    
    # ******* condition 5 ***********
    # EIF line 5
    
    b.con5_funs$u <- true_funcs$u
    b.con5_funs$B1 <- true_funcs$B1
    b.con5_funs$B2 <- true_funcs$B2
    b.con5_funs$B3 <- true_funcs$B3
    
    b.con5_funs$prob1.T1 <- true_funcs$prob1.T1; b.con5_funs$prob0.T1 <- true_funcs$prob0.T1
    
    # ******* condition 6 ***********
    
    # EIF line 6
    
    b.con6_funs$u <- true_funcs$u
    b.con6_funs$B1 <- true_funcs$B1
    b.con6_funs$B2 <- true_funcs$B2
    b.con6_funs$B3 <- true_funcs$B3
    b.con6_funs$B4 <- true_funcs$B4
    
    # ******* output ***********
    out <- list(b.con1_funs = b.con1_funs, b.con2_funs = b.con2_funs, b.con3_funs = b.con3_funs,
                b.con4_funs = b.con4_funs, b.con5_funs = b.con5_funs, b.con6_funs = b.con6_funs)
    
  }
  
  
  return(out)
}
  
  
  
component.T1T2 <- function(data, estimation_funs_list){
  
  list2env(estimation_funs_list, envir = environment())
  
  indicator.list <- list(
    c((data$T1 == 0)&(data$T2 == 0)),
    c((data$T1 == 1)&(data$T2 == 1)),
    c((data$T1 == 0)&(data$T2 == 0)),
    c((data$T1 == 1)),
    c((data$T1 == 0)),
    1
  )
  
  
  prob.list = list(prob0.T1*prob0.T2,
                   prob1.T1*prob1.T2,
                   prob0.T1*prob0.T2,
                   prob1.T1,
                   prob0.T1,
                   1)
  
  
  density.ratio.list = list(density.ratio.L1*density.ratio.L2,
                            (1/density.ratio.M1)*(1/density.ratio.M2),
                            density.ratio.L1,
                            (1/density.ratio.M1),
                            1,
                            1
  )
  
  
  # OR formula for joint distribution
  bayes.ratio.list = list(bayes.ratio.T1_L1 * (1/bayes.ratio.T1_M1)*
                            bayes.ratio.T2_L2 * bayes.ratio.T1_L2*
                            (1/bayes.ratio.T2_M2) * (1/bayes.ratio.T1_M2),
                          
                          # row 2
                          (1/bayes.ratio.T1_M1) * (prob1.T1/prob0.T1)*
                            (1/bayes.ratio.T2_M2) * (1/bayes.ratio.T1_M2)*
                            (prob1.T2/prob0.T2) * bayes.ratio.T1_L1,
                          
                          # row 3
                          bayes.ratio.T1_L1 * (1/bayes.ratio.T1_M1),
                          # row 4
                          (1/bayes.ratio.T1_M1) * (prob1.T1/prob0.T1),
                          
                          # row 5
                          1,
                          
                          # row 6
                          1
  )
  
  seq_diff = list(Y-u, u-B1, B1-B2, B2-B3, B3-B4, B4)
  
  return(
    list(indicator.list = indicator.list,
         prob.list = prob.list,
         density.ratio.list = density.ratio.list,
         bayes.ratio.list = bayes.ratio.list,
         seq_diff = seq_diff)
  )
}



longT_one_step <- function(component.list, type = "bayes"){
  list2env(component.list, envir = environment())
  
  if(type == "bayes"){ratio.list = bayes.ratio.list}else{ratio.list  = density.ratio.list}
  
  partial_results<-suppressMessages(pmap_dfc(list(indicator.list, prob.list, ratio.list,  seq_diff), 
                            function(indicator, prob, ratio,  seq_diff){
                              indicator/prob*ratio*seq_diff}
                            ) )
  
  counterfactual_outcome = rowSums(partial_results, na.rm = T)
  
  out <- data.frame(value = mean(counterfactual_outcome),
                    SE = sqrt(var(counterfactual_outcome)/length(counterfactual_outcome)),
                    EIF_var = var(counterfactual_outcome),
                    type = type
                    )
  
  return(list(out = out, counterfactual_outcome = counterfactual_outcome))
}



