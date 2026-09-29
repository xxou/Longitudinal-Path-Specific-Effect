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
# 1 10101010;101010;1010;10 vs 00000000;000000;0000;00 -0.6531351 0.2730128 -1.1882305 -0.1180398  0.0167
# 2 11111111;111111;1111;11 vs 00000000;000000;0000;00  0.1265960 0.5288754 -0.9099808  1.1631728  0.8108
# 3 11111111;111111;1111;11 vs 10101010;101010;1010;10  0.7797311 0.4703778 -0.1421925  1.7016547  0.0974


# GLM
# active    Effect        SE   CI.lower   CI.upper P.value
# 1 10101010;101010;1010;10 vs 00000000;000000;0000;00 -1.824171 0.6906523  -3.177825 -0.4705175  0.0083
# 2 11111111;111111;1111;11 vs 00000000;000000;0000;00 -4.654182 6.0601860 -16.531928  7.2235648  0.4425
# 3 11111111;111111;1111;11 vs 10101010;101010;1010;10 -2.830010 5.6780294 -13.958743  8.2987228  0.6182
