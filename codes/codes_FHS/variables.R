# ************************************
# summarize variables
# ************************************

# Baseline covariates
baseline_cov <- c("SEX", "AGE1", "HGT1", "EDUC", "OCCUP", "MARITAL", "TOBACCO")

# Time varying L
L_exam2 <- c("SBP1", "TC2")
L_exam4 <- c("SBP2", "TC4", "HRX4")
L_exam5 <- c("SBP4", "TC5", "HRX5")
L_exam7 <- c("SBP5", "TC7", "HRX7")

# Mediator M
M_exam2 <- "dBMI2"
M_exam4 <- "dBMI4"
M_exam5 <- "dBMI5"
M_exam7 <- "dBMI7"

# Exposure A
A_exam2 <- "SMOKE2"
A_exam4 <- "SMOKE4"
A_exam5 <- "SMOKE5"
A_exam7 <- "SMOKE7"

# outcome
Y_outcome = 'SBP7'

