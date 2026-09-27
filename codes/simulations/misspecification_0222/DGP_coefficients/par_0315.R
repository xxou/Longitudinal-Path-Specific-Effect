# ***************************************
# longitudinal treatments
# collaborate with Xinwei
# Mar - 15 - 2025
# Simplest data for misspecification
# ***************************************

# paste(round(runif(n=5, -0.5,0.5),2), collapse = ", ") T1
# paste(round(runif(n=8, -0.5,0.5),2), collapse = ", ") T2
# paste(round(runif(10,-0.2,0.5),2), collapse = ", ")  M_i
# paste(round(runif(11,-0.5,1),2), collapse = ", ")   Y

# high correlation between T1 and T2

par_T1 <- c(-0.21, -0.25, 0.35, -0.26, 0.49)
par_M1 <- c(0.5, -0.02, 0.32, 0.13, -0.14, 0.15)
par_L1 <- c(0.33, -0.14, 0.21, 0.14, 0.46, 0.14, 0.34)
par_T2 <- c(-0.04, 0.13, 0.26, -0.16, -0.38, 0.41, 0.12, -0.29)
par_M2 <- c(0.24, -0.06, 0.29, 0.49, -0.15, 0.3, -0.13, 0.18, 0.48)
par_L2 <- c(-0.11, 0.01, 0.46, 0.12, -0.04, -0.06, 0.24, 0.43, 0.45, 0.18)
par_Y <- c(0.19, 0.55, 0.6, 0.2, 0.14, 0.65, 0.15, -0.37, 0.54, 0.66, 0.9)

# 2.204582; 13.08338
# value          SE      EIF_var    type
# 1 2.199669 0.001142514 13.05339   bayes
# 2 2.199668 0.001142518 13.05348 density