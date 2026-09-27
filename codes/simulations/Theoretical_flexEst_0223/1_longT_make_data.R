
# ***************************************
# longitudinal treatments
# collaborate with Xinwei
# Feb - 23 - 2025
# ***************************************

# ***************************************
#            1. make data   -------------

longT_make_data <- function(n=2000, par_T1, par_M1, par_L1, par_T2, par_M2, par_L2, par_Y){

  X1 <- rbinom(n, size = 1, prob = 0.5)
  X2 <- runif(n, -1, 1)
  X3 <- rnorm(n, mean = 0, sd = 1)
  
  # T1
  # fmla_T1_intr <- T1~X1+X2+X3+I(X2^2)+X2:X3; hist(T1); hist(prob); summary(prob)
  T1 = rbinom(n, 1, prob = plogis(cbind(X1, X2, X3, X2^2, X2*X3) %*% par_T1))
  
  t1_m1 = t1_l1 = t1_m2 = t1_l2 = t1_y =  T1
  
  # M1 
  # hist(M1);
  M1 = cbind(X1, X2, X3, t1_m1, X2^2, X1*X3, t1_m1*X3) %*% par_M1 + rnorm(n, 0, 1)

  # L1 
  # hist(L1);
  L1 = cbind(X1, X2, X3, t1_l1, M1, X2^2) %*% par_L1 + rnorm(n, 0, 1)
  
  
  # T2 
  # hist(T2);  hist(prob); summary(prob)
  T2 = rbinom(n, 1, prob = plogis(cbind(X1, X2, X3, T1, M1, L1, X3^2, X1*X3) %*% par_T2)) 
  
  t2_m2 = t2_l2 = t2_y = T2
  
  # M2 
  # hist(M2);
  # par_M2 <- c(-0.6, 0.2, -0.5, -1.0, -0.7, 0.2, -0.7, 0.2) # changed
  M2 = cbind(X1, X2, X3, t1_m2, M1, L1, t2_m2, X2^2) %*% par_M2 + rnorm(n, 0, 1)
  
  # L2 
  # hist(L2);
  L2 = cbind(X1, X2, X3, t1_l2, M1, L1, t2_l2, M2, X3^2, X2*X3) %*% par_L2 + rnorm(n, 0, 1)
  
  # Y 
  # hist(Y)
  Y = cbind(X1, X2, X3, t1_y, M1, L1, t2_y, M2, L2, X3^2) %*% par_Y + rnorm(n, 0, 1)
  
  data = data.frame(X1, X2, X3, T1, M1, L1,T2, M2, L2, Y)
  
  return(data)
}



# ***************************************
#            2. get truth   -------------

get_truth_longT <- function(n = 1000, T1 = c(), T2 = c(), par_T1, par_M1, par_L1, par_T2, par_M2, par_L2, par_Y){
  # covariates
  X1 <- rbinom(n, size = 1, prob = 0.5)
  X2 <- runif(n, -1, 1)
  X3 <- rnorm(n, mean = 0, sd = 1)
  
  t1_m1 = T1[1]; t1_l1 = T1[2]; t1_m2 = T1[3]; t1_l2 = T1[4]; t1_y = T1[5]
                                t2_m2 = T2[1]; t2_l2 = T2[2]; t2_y = T2[3]
  
  
  # M1 
  # hist(M1); hist(prob)
  M1 = cbind(X1, X2, X3, t1_m1, X2^2, X1*X3, t1_m1*X3) %*% par_M1 + rnorm(n, 0, 1)
  # L1 
  L1 = cbind(X1, X2, X3, t1_l1, M1, X2^2) %*% par_L1 + rnorm(n, 0, 1)
  # M2 
  # par_M2 <- c(-0.6, 0.2, -0.5, -1.0, -0.7, 0.2, -0.7, 0.2) # changed
  M2 = cbind(X1, X2, X3, t1_m2, M1, L1, t2_m2, X2^2) %*% par_M2 + rnorm(n, 0, 1)
  # L2 
  L2 = cbind(X1, X2, X3, t1_l2, M1, L1, t2_l2, M2, X3^2, X2*X3) %*% par_L2 + rnorm(n, 0, 1)
  # Y 
  Y = cbind(X1, X2, X3, t1_y, M1, L1, t2_y, M2, L2, X3^2) %*% par_Y + rnorm(n, 0, 1)
  
  out = data.frame(mean = mean(Y), var = var(Y))
  
  return(out)
}

# get_truth_longT(1e7, T1 = c(0, 1, 0, 1, 0), T2 = c(0, 1, 0))
#    mean     var
# 1 1.590974 10.11683



