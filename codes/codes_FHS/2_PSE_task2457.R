# ***********************************************
# PSE for FHS data
# {flexPaths} package
# ***********************************************

freshr::freshr()

# load package
pacman::p_load(flexPaths, dbarts, SuperLearner, xgboost, ranger, glmnet, kernlab)

source(here::here("codes/functions/SL.learner.R"))
source(here::here("codes/codes_FHS/variables.R"))

# load data
load(here::here("data/FHS_data_use2457.RData"))

data = FHS_data_use
# ********************************
# PSE run            -------------
SL_library = c('SL.earth','SL.xgboost','SL.ranger.1000', 'SL.dbarts', 'SL.glm', 'SL.ksvm' ) 
# SL_library = c('SL.glm')


FHS_fit<- pathsInfo(data = data,
                    Y = Y_outcome,
                    A = c(A_exam2, A_exam4, A_exam5, A_exam7),
                    cov_x = c(baseline_cov, L_exam2),
                    M.list = list(
                     M_exam2, L_exam4,
                     M_exam4, L_exam5, 
                     M_exam5, L_exam7,
                     M_exam7
                    ),
                    estimation = "EIF",
                    model.treatment = ~ SuperLearner(SL.library = SL_library, family = binomial()),
                    model.outcome =  ~ SuperLearner(SL.library = SL_library, family = gaussian()),
                    model.iter = ~ SuperLearner(SL.library = SL_library, family = gaussian())
)

mp_m<-mflexPotential(active = list(a1=c(1,0,1,0,1,0,1,0),
                                  a2=c(NA,NA,1,0,1,0,1,0),
                                  a3=c(NA,NA,NA,NA,1,0,1,0),
                                  a4=c(NA,NA,NA,NA,NA,NA,1,0)),
                    FHS_fit)
mp1<-mflexPotential(active = list(a1=c(1,1,1,1,1,1,1,1),
                                  a2=c(NA,NA,1,1,1,1,1,1),
                                  a3=c(NA,NA,NA,NA,1,1,1,1),
                                  a4=c(NA,NA,NA,NA,NA,NA,1,1)),
                    FHS_fit)
mp0<-mflexPotential(active = list(a1=c(0,0,0,0,0,0,0,0),
                                  a2=c(NA,NA,0,0,0,0,0,0),
                                  a3=c(NA,NA,NA,NA,0,0,0,0),
                                  a4=c(NA,NA,NA,NA,NA,NA,0,0)),
                    FHS_fit)

Effect = flexEffect(p1 = list(mp1, mp1, mp_m), p0 = list(mp0, mp_m, mp0), scale = "diff")
Effect

save(Effect, file = here::here("output/hpc_data/FHS_PSE_Effect2457.RData"))


load(here::here("output/hpc_data/FHS_PSE_Effect2457.RData"))

# SL
# active     Effect        SE   CI.lower   CI.upper P.value
# 1 10101010;101010;1010;10 vs 00000000;000000;0000;00 -0.8399249 0.3116063 -1.4506621 -0.2291878  0.0070
# 2 11111111;111111;1111;11 vs 00000000;000000;0000;00  0.1477803 0.5872147 -1.0031393  1.2987000  0.8013
# 3 11111111;111111;1111;11 vs 10101010;101010;1010;10  0.9877053 0.7158857 -0.4154049  2.3908155  0.1677


# GLM
# active    Effect        SE   CI.lower   CI.upper P.value
# 1 10101010;101010;1010;10 vs 00000000;000000;0000;00 -1.983395 0.7840341  -3.520073 -0.4467161  0.0114
# 2 11111111;111111;1111;11 vs 00000000;000000;0000;00 -5.354443 6.5343891 -18.161610  7.4527243  0.4125
# 3 11111111;111111;1111;11 vs 10101010;101010;1010;10 -3.371048 6.0117231 -15.153809  8.4117125  0.5750
