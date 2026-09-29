# ************************************************
# clean data from FHS data
# exam 2, 4, 5, 7
# 2024-Jun-24
# ************************************************

freshr::freshr()


# load package
pacman::p_load(tidyverse, gtsummary)


# dataset
ex0_7d <- read.csv("data/FHS_Data/Datasets/CSV/ex0_7d.csv")
vr_wkthru_ex32_0_0997d <- read.csv("data/FHS_Data/Datasets/CSV/vr_wkthru_ex32_0_0997d.csv")
vr_survdth_2018_a_1268d <- read.csv("data/FHS_Data/Datasets/CSV/vr_survdth_2018_a_1268d.csv")
vr_survcvd_2018_a_1267d <- read.csv("data/FHS_Data/Datasets/CSV/vr_survcvd_2018_a_1267d.csv")
vr_cancer_2019_a_1162d <- read.csv("data/FHS_Data/Datasets/CSV/vr_cancer_2019_a_1162d.csv")


# ********************************************
#  1.coding variables                  -------
# ********************************************
# excel: ou_codebook_clean
# rename variables
# modify levels


## (1)  select variables   ----------
use_ex0_7d <- ex0_7d |> 
  select(PID, MF5, MF447, MF79, MF71, MF522)
use_vr_wkthru <- vr_wkthru_ex32_0_0997d |> 
  select(PID,
         SEX, 
         matches("^AGE[1-7]$"),
         HGT1, HGT4, HGT5,
         matches("^SBP[1-7]$"),
         matches("^TC[1-7]$"),
         matches("^HRX[1-7]$"),
         matches("^CPD[1-7]$"),
         matches("^CURRSMK[1-7]$"),
         matches("^BMI[1-7]$"),
         matches("^WGT[1-7]$"),
         matches("^DATE[1-7]$"),
         matches("^ATT[1-7]$"),
         DMRX1
  )
use_vr_survdth <- vr_survdth_2018_a_1268d |> 
  select(PID, DATEDTH)
use_vr_survcvd <- vr_survcvd_2018_a_1267d |> 
  select(PID, CVD, CVDDATE)
use_vr_cancer <- vr_cancer_2019_a_1162d |> 
  select(PID, D_DATE)


# FHS data
FHS_data_basic <- use_ex0_7d |> 
  full_join(use_vr_wkthru, by = join_by(PID)) |> 
  full_join(use_vr_survdth, by = join_by(PID)) |> 
  full_join(use_vr_survcvd, by = join_by(PID)) 



## (2)  variables  labels and levels  ----------
FHS_data0 <- FHS_data_basic |> 
  mutate( 
    ####  gender
    SEX = factor(SEX, labels = c("Male", "Female")),
    
    #### BMI
    HGT2 = HGT1, HGT3 = HGT1, 
    HGT6 = HGT5, HGT7 = HGT5,
    across(matches("^HGT[1-7]$"), ~ .x * 0.0254, .names = "{.col}_m"),
    across(matches("^WGT[1-7]$"), ~ .x * 0.4535, .names = "{.col}_kg"),
    
    dBMI1 = WGT1_kg / (HGT1_m)^2,
    dBMI2 = WGT2_kg / (HGT2_m)^2,
    dBMI3 = WGT3_kg / (HGT3_m)^2,
    dBMI4 = WGT4_kg / (HGT4_m)^2,
    dBMI5 = WGT5_kg / (HGT5_m)^2,
    dBMI6 = WGT6_kg / (HGT6_m)^2,
    dBMI7 = WGT7_kg / (HGT7_m)^2,
    
    
    #### education
    EDUC = case_when(MF5 <= 3 ~ "< High school",
                     MF5 %in% c(4, 5) ~ "High school",
                     MF5 %in% c(6, 7, 8, 10) ~ "college or higher"
    ),
    EDUC = factor(EDUC, levels = c("< High school", "High school", "college or higher")),
    
    
    #### occupation
    OCCUP = case_when(MF447 == 1 ~ "Professional",
                      (MF447 == 2) | (MF447 == 3) ~ "Executive/Supervisory",
                      MF447 == 4 ~ "Technical",
                      MF447 == 5 ~ "Laborer",
                      MF447 == 6 ~ "Clerical",
                      MF447 == 7 ~ "Sales",
                      MF447 == 8 ~ "Housewife"
    ),
    
    
    #### marital status
    MARITAL = case_when(MF79 == 1 ~ "Single",
                        MF79 <= 4 ~ "Married",
                        MF79 <= 6 ~ "divorce or widowed"),
    MARITAL = factor(MARITAL, levels = c("Single", "Married", "divorce or widowed")),
    
    #### tobacco 
    TOBACCO = case_when(MF71 %in% c(5, 7) ~ "Nonsmokers",
                        MF71 %in% c(1, 3) ~ "Smokers",
                        MF71 == 2 ~ "Quitters"
    ),
    TOBACCO = factor(TOBACCO, levels = c("Quitters", "Nonsmokers", "Smokers")),
    
    #### SMOKE 
    across(starts_with("CPD"), ~ if_else(.x > 0, 1, 0), .names = "SMOKE{sub('CPD', '', .col)}")
    # no SMOKE6
  ) |> rename(CURR_DIABETES1 = "DMRX1")


FHS_data_all <- FHS_data0 |> select(
  PID, 
  # baseline
  SEX, AGE1, HGT1, EDUC, OCCUP, MARITAL, TOBACCO, 
  # covariates L and outcome Y
  contains('SBP'), contains('TC'), contains('HRX'),
  # exposure
  contains('SMOKE'),
  # mediator
  contains('dBMI'),
  # validation
  matches("^DATE[1-7]$"), matches("^ATT[1-7]$"), DATEDTH,
  CURR_DIABETES1, CVDDATE, CVD
)



## (3)  update missing information  ----------
tbl_summary(FHS_data_all, include = - PID,
            statistic = list(all_continuous() ~ '{median} ({min}, {max})',
                             all_categorical() ~ "{n} ({p}%)"),
            digits = list(all_continuous() ~ c(2, 2),
                          all_categorical() ~ c(0, 1))
)

# For missing BMi value or smoking status at a single time period, 
# we carry forward the last observed value/status for one exam period only.

FHS_data_all_update <- FHS_data_all |> mutate(
  dBMI7 = coalesce(dBMI7, dBMI6),
  dBMI6 = coalesce(dBMI6, dBMI5),
  dBMI5 = coalesce(dBMI5, dBMI4),
  dBMI4 = coalesce(dBMI4, dBMI3),
  dBMI3 = coalesce(dBMI3, dBMI2),
  dBMI2 = coalesce(dBMI2, dBMI1),
  
  # no SMOKE6
  SMOKE7 = coalesce(SMOKE7, SMOKE5),
  SMOKE5 = coalesce(SMOKE5, SMOKE4),
  SMOKE4 = coalesce(SMOKE4, SMOKE3),
  SMOKE3 = coalesce(SMOKE3, SMOKE2),
  SMOKE2 = coalesce(SMOKE2, SMOKE1)
)

tbl_summary(FHS_data_all_update, include = - PID)
# timeline
# baseline [exam 1]
# exam 2 --> exam 4 --> exam 5 --> exam 7


# ********************************************
#  2.  exclusion                       -------
# ********************************************
FHS_data_all_update |> names()
mean_DATE2 <- mean(FHS_data_all_update$DATE2, na.rm = T)

cancer_baseline <- use_vr_cancer |> 
  left_join( FHS_data_all_update[ , c('PID', 'DATE2')], by = "PID") |> 
  mutate(DATE2 = if_else(is.na(DATE2), mean_DATE2, DATE2)) |> 
  filter(D_DATE < DATE2) |> distinct(PID)


FHS_data_inclusion <- FHS_data_all_update |> 
  # (1) death or loss to followup during the period before exam 7
  # exam 1-7
  mutate(
    DEATH_ex = if_else(
      !is.na(DATEDTH) & DATEDTH < DATE7, 1, 0),
    ATT_ex = if_else(if_all(c(ATT1, ATT2, ATT4, ATT5, ATT7), ~ .x == 1), 1, 0) ) |> 
  filter(DEATH_ex == 0, ATT_ex == 1) |> 
  
  # (2) diagnosis of diabetes, cancer, or CVD at baseline
  filter(!(PID %in% cancer_baseline$PID),
         CURR_DIABETES1 == 0,
         !(CVDDATE < DATE2) ) |> 
  
  select(-c(PID, SBP3, SBP6, 
            TC1, TC3, TC6, 
            HRX6,
            SMOKE1, SMOKE3,
            dBMI1, dBMI3, dBMI6, 
            contains("DATE"), contains("ATT"),
            CURR_DIABETES1, CVD, DEATH_ex))  


# (3) missing data
FHS_data_use <- FHS_data_inclusion |>  drop_na()


tbl_summary(FHS_data_use)

tbl_summary(FHS_data_use, by = TOBACCO) |> 
  add_overall()|> 
  add_p() 

# ********************************************
#      save data                       -------
# ********************************************
save(FHS_data_use, file = "data/FHS_data_use2457.RData")
save(FHS_data_all, FHS_data_all_update, FHS_data_inclusion,  FHS_data_use, file = "data/FHS_data2457.RData")
