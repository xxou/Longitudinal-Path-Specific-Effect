
# ***************************************
# longitudinal treatments
# collaborate with Xinwei
# Feb - 22 - 2025
# Simplest data for misspecification
# ***************************************


# ***************************************
#            1. make data   -------------

longT_make_simplest_data <- function(n=2000, par_T1, par_M1, par_L1, par_T2, par_M2, par_L2, par_Y){

  X <- matrix(c(runif(n,0,1), runif(n, 0,1),runif(n, 0,1),runif(n, 0,1)), byrow = F, ncol = 4)
  
  
  
  # T1
  # hist(prob); summary(prob); hist(T1)
  T1 = rbinom(n, 1, prob = plogis(cbind(1, X) %*% par_T1))
  
  t1_m1 = t1_l1 = t1_m2 = t1_l2 = t1_y =  T1
  
  # M1 
  # hist(M1); 
  M1 = cbind(1, X, t1_m1) %*% par_M1 + rnorm(n, 0, 1)

  # L1 
  # hist(L1)
  L1 = cbind(1, X, t1_l1, M1) %*% par_L1 + rnorm(n, 0, 1)
  
  
  # T2 
  # hist(T2); hist(prob); summary(prob)
  T2 = rbinom(n, 1, prob = plogis(cbind(1, X, T1, M1, L1) %*% par_T2)) 
  
  t2_m2 = t2_l2 = t2_y = T2
  
  # M2 
  # hist(M2)
  M2 = cbind(1, X, t1_m2, M1, L1, t2_m2) %*% par_M2 + rnorm(n, 0, 1)
  
  # L2
  # hist(L2)
  L2 = cbind(1, X, t1_l2, M1, L1, t2_l2, M2) %*% par_L2 + rnorm(n, 0, 1)
  
  # Y 
  # hist(Y)
  Y = cbind(1, X, t1_y, M1, L1, t2_y, M2, L2) %*% par_Y + rnorm(n, 0, 1)
  
  data = data.frame(X, T1, M1, L1, T2, M2, L2, Y)
  
  return(data)
}



# ***************************************
#            2. get truth   -------------

get_truth_longT <- function(n = 1000, T1 = c(), T2 = c(), par_T1, par_M1, par_L1, par_T2, par_M2, par_L2, par_Y){
  # covariates
  X <- matrix(c(runif(n,0,1), runif(n, 0,1),runif(n, 0,1),runif(n, 0,1)), byrow = F, ncol = 4)
  
  t1_m1 = T1[1]; t1_l1 = T1[2]; t1_m2 = T1[3]; t1_l2 = T1[4]; t1_y = T1[5]
  t2_m2 = T2[1]; t2_l2 = T2[2]; t2_y = T2[3]
  
  
  # M1 
  # hist(M1); 
  M1 = cbind(1, X, t1_m1) %*% par_M1 + rnorm(n, 0, 1)
  
  # L1 
  # hist(L1)
  L1 = cbind(1, X, t1_l1, M1) %*% par_L1 + rnorm(n, 0, 1)
  
  # M2 
  # hist(M2)
  M2 = cbind(1, X, t1_m2, M1, L1, t2_m2) %*% par_M2 + rnorm(n, 0, 1)
  
  # L2
  # hist(L2)
  L2 = cbind(1, X, t1_l2, M1, L1, t2_l2, M2) %*% par_L2 + rnorm(n, 0, 1)
  
  # Y 
  # hist(Y)
  Y = cbind(1, X, t1_y, M1, L1, t2_y, M2, L2) %*% par_Y + rnorm(n, 0, 1)
  
  out = data.frame(mean = mean(Y), var = var(Y))
  
  return(out)
}

# get_truth_longT(1e8, T1 = c(0, 1, 0, 1, 0), T2 = c(0, 1, 0), 
#                 par_T1, par_M1, par_L1, par_T2, par_M2, par_L2, par_Y)


