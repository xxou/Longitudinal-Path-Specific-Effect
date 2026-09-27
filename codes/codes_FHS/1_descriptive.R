# ***********************************************
# Descriptive analysis for FHS data
# ***********************************************

freshr::freshr()

# load package
pacman::p_load(gtsummary, gt, tidyverse)

source(here::here("codes/codes_FHS/variables.R"))

# load data
load(here::here("data/FHS_data2457.RData"))



# ***************************
# check variables   ---------

tbl_summary(FHS_data_all, include = - PID,
            statistic = list(all_continuous() ~ '{median} ({min}, {max})',
                             all_categorical() ~ "{n} ({p}%)"),
            digits = list(all_continuous() ~ c(2, 2),
                          all_categorical() ~ c(0, 1))
)

tbl_summary(FHS_data_use, by = TOBACCO,
            
            statistic = list(all_continuous() ~ '{mean} ± {sd})',
                             all_categorical() ~ "{n} ({p}%)"),
            digits = list(all_continuous() ~ c(1, 1),
                          all_categorical() ~ c(0, 1))
            ) |> 
  add_overall() 


FHS_data_use |>  
  tbl_summary(by = TOBACCO,
            include = c(SEX, AGE1, HGT1, EDUC, OCCUP, MARITAL, SBP1), 
            label = list(
              SEX ~ "Sex",
              AGE1 ~ "Age (years)",
              HGT1 ~ "Height (m)",
              EDUC ~ "Education level",
              OCCUP ~ "Occupation",
              MARITAL ~ "Marital status",
              SBP1 ~ "SBP (mm Hg)"
            ),
            statistic = list(all_continuous() ~ '{mean} ± {sd}',
                             all_categorical() ~ "{n} ({p}%)"),
            digits = list(all_continuous() ~ c(1, 1),
                          all_categorical() ~ c(0, 1))
) |> 
  add_overall()  -> tab1
tab1 |> as_gt() |> as_latex() -> tab1_gt
cat(tab1_gt)
