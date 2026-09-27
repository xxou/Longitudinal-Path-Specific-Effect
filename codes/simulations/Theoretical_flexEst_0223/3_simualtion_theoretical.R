# *****************************************************************
# local parallel: One step estimator for longitudinal treatments
# *****************************************************************
# DGP:  Feb - 23 - 2025
# run locally
# Xiaxian Ou
# two treatments T1 -> M1 -> L1 -> T2 -> M2 -> L2 -> Y
# *****************************************************************

# *****************************************************************
#                       1 simulation                 --------------
# *****************************************************************

freshr::freshr()

# package
library(tidyverse)

# run process
source("codes/simulations/Theoretical_flexEst_0223/1_longT_make_data.R")
source("codes/functions/est_funs_GLM.R")
source("codes/functions/one_step_estimator.R")


# ***********************************
## one simulation      --------------

run_simulation <- function(seed, n){
  
  # set seed
  set.seed(seed)
  cat("seed:", seed, "; n:",  n, "\n")
  
  source("codes/simulations/Theoretical_flexEst_0223/DGP_coefficients/par_0323.R")  
  data = longT_make_data(n, par_T1, par_M1, par_L1, par_T2, par_M2, par_L2, par_Y)
  
  # ******** correct specification **************
  # formula
  f.pro.T.list = list(pro.T1_formula = T1~X1+X2+X3+I(X2^2)+X2*X3, 
                      pro.T2_formula = T2~X1+X2+X3+T1+M1+L1+I(X3^2)+X1*X3)
  
  f.mediator.list = list(M1_formula = M1~X1+X2+X3+T1+I(X2^2)+X1:X3+T1:X3,
                         L1_formula = L1~X1+X2+X3+T1+M1+I(X2^2),
                         M2_formula = M2~X1+X2+X3+T1+M1+L1+T2+I(X2^2),
                         L2_formula = L2~X1+X2+X3+T1+M1+L1+T2+M2+I(X3^2)+X2:X3)
  
  
  f.bayes.list = list(T1_M1_formula = T1~X1+X2+X3+I(X2^2)+I(X3^2)+X1*X3+X2*X3+I(X2^2)*X3+X1*I(X3^2)+M1+M1*X3,
                      T1_L1_formula = T1~X1+X2+X3+I(X2^2)+I(X3^2)+X1*X3+X2*X3+I(X2^2)*X3+X1*I(X3^2)+M1+M1*X3+L1,
                      
                      T1_M2_formula = T1~X1+X2+X3+I(X2^2)+I(X3^2)+X1*X3+X2*X3+I(X2^2)*X3+X1*I(X3^2)+M1+M1*X3+L1+M2+T2,
                      T1_L2_formula = T1~X1+X2+X3+I(X2^2)+I(X3^2)+X1*X3+X2*X3+I(X2^2)*X3+X1*I(X3^2)+M1+M1*X3+L1+M2+T2+L2,
                      
                      T2_M2_formula = T2~X1+X2+X3+T1+M1+L1+I(X3^2)+X1*X3+M2+I(X2^2),
                      T2_L2_formula = T2~X1+X2+X3+T1+M1+L1+I(X3^2)+X1*X3+M2+I(X2^2)+L2+X2*X3)
  
  
  seq.reg.list <- list(
    # E(Y| ....):
    u = predY ~ X1+X2+X3+M1+L1+M2+L2+I(X3^2)+T1+T2,
    # B1; consider E(L2|...) in addition:
    B1 = predY ~ X1+X2+X3+M1+L1+M2+I(X3^2)+X2*X3+T1+T2,
    # B2; consider E(M2|...) in addition:
    B2 = predY ~ X1+X2+X3+M1+L1+I(X3^2)+X2*X3+I(X2^2)+T1+T2,
    # B3; consider E(L1|...) in addition:
    B3 = predY ~ X1+X2+X3+M1+I(X3^2)+X2*X3+I(X2^2)+T1,
    # B4; consider E(M1|...) in addition
    B4 = predY ~ X1+X2+X3+I(X3^2)+X2*X3+I(X2^2)+X1*X3+T1*X3+T1
  )  
  
  
  funcs <- GLM_nuisance_funs(data, f.pro.T.list, f.mediator.list, f.bayes.list, seq.reg.list, joint = "OR", At_correction = TRUE)
  
  # ******** results **************
  component <- component.T1T2(data, funcs)
  estimate = rbind(longT_one_step(component, type = "density")$out, 
                   longT_one_step(component, type = "bayes")$out)
  estimate$seed = seed
  estimate$sample_size = n
  
  return(estimate)
}

# **************************************
## parallel simulation -----------------

library(furrr)

plan(multisession)
# setting  ********************
unit_n = 1000

param_grid <- expand.grid(seed = 1:unit_n,
                          n = c(500, 1000, 2000, 4000, 8000)
)
# 
# results <- future_map2(param_grid$seed, param_grid$n, run_simulation, 
#                        .options = furrr_options(seed = TRUE))

nbreaks = 10
results = vector('list', nbreaks)
param_grid_splits <- split(param_grid, cut(seq_len(nrow(param_grid)), breaks = nbreaks, labels = FALSE))

for (i in seq_len(nbreaks)) {
  
  cat("i:", i, " ------------------------ ", "\n")
  
  results[[i]] <- future_map2(param_grid_splits[[i]]$seed, param_grid_splits[[i]]$n, 
                              run_simulation, .options = furrr_options(seed = NULL))
}

# Combine results into a dataframe
results <- bind_rows(results)
save(results, file =  "output/local_parallel/correct/correct_specify0323_At.Rdata")


# *****************************************************************
#                       2 summarize                 --------------
# *****************************************************************

freshr::freshr()
pacman::p_load(tidyverse,gtsummary,gt, ggpubr)

load("output/local_parallel/correct/correct_specify0323_At.Rdata")

simu.longT <- results
head(simu.longT)
# summary
simu.longT.summary <- simu.longT |> 
  mutate( true_value = 0.394) |> 
  mutate(
    error = value - true_value,
    square_error = (value - true_value)^2,
    
    CI.upper = value + 1.96*SE,
    CI.lower = value - 1.96*SE,
    
    cover = if_else( (true_value <= CI.upper) & (true_value >= CI.lower), 1, 0),
    length = CI.upper - CI.lower)

simu.longT.evaluate <- simu.longT.summary |>   
  group_by(sample_size, type) |>  
  mutate(dev = value - mean(value),
         square_dev = (value - mean(value))^2) |> 
  summarise(
    bias=mean(error),
    se_bias = sqrt(mean(square_dev)/n()),
    mse = mean(square_error),
    
    sd_est = sd(value),
    var_est = mean(square_dev),
    se_var_est = sqrt(var(square_dev)/n()),
    
    cover_rate = mean(cover,na.rm=T),
    avg_length =  mean(length,na.rm=T)
    
  ) |>  ungroup() |>  
  mutate(
    scaled_bias = sqrt(sample_size)*bias,
    se_scaled_bias = sqrt(sample_size)*se_bias,
    CI.up.scaled_bias = scaled_bias+1.96*se_scaled_bias,
    CI.low.scaled_bias = scaled_bias-1.96*se_scaled_bias,
    
    scaled_var_est = sample_size*var_est,
    se_scaled_var_est = sample_size*se_var_est,
    CI.up.scaled_var_est = scaled_var_est+1.96*se_scaled_var_est,
    CI.low.scaled_var_est = scaled_var_est-1.96*se_scaled_var_est,
    
    true_var = 13.94
  )



### (1) display table 
gt.simu.longT.evaluate <- simu.longT.evaluate |> 
  select(type, sample_size, bias, sd_est, mse, cover_rate, avg_length) |> 
  # pivot_wider(
  #   names_from = model,
  #   values_from = c(bias, sd_est, mse, cover_rate, avg_length),
  #   names_glue = "{.value}_{model}" 
  # ) |> 
  gt(groupname_col = "type", row_group_as_column = F) |> 
  tab_spanner(label = "Bias", columns = matches("^bias_"))|> 
  tab_spanner(label = "Standard Error", columns = matches("^sd_est_")) |> 
  tab_spanner( label = "Mean Squared Error",columns = matches("^mse_")) |> 
  tab_spanner(label = "Coverage Rate", columns = matches("^cover_rate_")) |> 
  tab_spanner(label = "Average Length",columns = matches("^avg_length_")) |> 
  fmt_number(columns = -sample_size, decimals = 3)
gt.simu.longT.evaluate

cat(as_latex(gt.simu.longT.evaluate))


### (2) plots 

simu.longT.evaluate.bayes <- simu.longT.evaluate |> filter(type == "bayes")
simu.longT.evaluate.density <- simu.longT.evaluate |> filter(type == "density")

ggplot(simu.longT.evaluate.bayes, aes(x = sample_size, y = scaled_bias)) +
  geom_point(color = "red", size = 2) +
  geom_line() +
  geom_errorbar(aes(ymin = CI.low.scaled_bias, ymax = CI.up.scaled_bias), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "blue") +
  labs(x = "Sample Size", y = expression(sqrt(n)-Bias)) +
  scale_y_continuous(limits = c(-0.5, 0.5))+
  theme_bw() +
  scale_x_continuous(
    breaks = unique(simu.longT.evaluate.bayes$sample_size),  # Set breaks at sample_size values
    labels = unique(simu.longT.evaluate.bayes$sample_size)  # Set labels as sample_size
  )->simu_rootn_bias_correct_bayes


ggplot(simu.longT.evaluate.density, aes(x = sample_size, y = scaled_bias)) +
  geom_point(color = "red", size = 2) +
  geom_line() +
  geom_errorbar(aes(ymin = CI.low.scaled_bias, ymax = CI.up.scaled_bias), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "blue") +
  labs(x = "Sample Size", y = expression(sqrt(n)-Bias)) +
  theme_bw() +
  scale_y_continuous(limits = c(-0.5, 0.5))+
  # facet_grid(rows = vars(model),scales = "fixed") +
  scale_x_continuous(
    breaks = unique(simu.longT.evaluate.density$sample_size),  # Set breaks at sample_size values
    labels = unique(simu.longT.evaluate.density$sample_size)  # Set labels as sample_size
  )->simu_rootn_bias_correct_density



# root n variance
ggplot(simu.longT.evaluate.bayes, aes(x = sample_size, y = scaled_var_est)) +
  geom_point(color = "red", size = 2) +
  geom_line() +
  geom_errorbar(aes(ymin = CI.low.scaled_var_est, ymax = CI.up.scaled_var_est), width = 0.2) +
  labs(x = "Sample Size", y = expression(n-Variance)) +
  theme_bw() +
  geom_hline(yintercept = simu.longT.evaluate.bayes$true_var, linetype = "dashed", color = "blue") +
  scale_y_continuous(limits = c(10, 17))+
  # facet_grid(rows = vars(model), scales = "fixed") +
  scale_x_continuous(
    breaks = unique(simu.longT.evaluate.bayes$sample_size),  # Set breaks at sample_size values
    labels = unique(simu.longT.evaluate.bayes$sample_size)  # Set labels as sample_size
  )->simu_n_variance_correct_bayes





ggplot(simu.longT.evaluate.density, aes(x = sample_size, y = scaled_var_est)) +
  geom_point(color = "red", size = 2) +
  geom_line() +
  geom_errorbar(aes(ymin = CI.low.scaled_var_est, ymax = CI.up.scaled_var_est), width = 0.2) +
  labs(x = "Sample Size", y = expression(n-Variance)) +
  theme_bw() +
  scale_y_continuous(limits = c(10, 17))+
  geom_hline(yintercept = simu.longT.evaluate.density$true_var, linetype = "dashed", color = "blue") +
  #facet_grid(rows = vars(model), scales = "fixed") +
  scale_x_continuous(
    breaks = unique(simu.longT.evaluate.density$sample_size),  # Set breaks at sample_size values
    labels = unique(simu.longT.evaluate.density$sample_size)  # Set labels as sample_size
  )->simu_n_variance_correct_density


ggarrange(simu_rootn_bias_correct_bayes, simu_n_variance_correct_bayes, nrow = 1) -> longT_correct_bayes
annotate_figure(
  longT_correct_bayes #,
 # top = text_grob("Bayes:At Correction", face = "bold", size = 14)
) ->  plot_longT_correct_bayes
plot_longT_correct_bayes

ggsave(plot_longT_correct_bayes, file="output/plots/plot_longT_correct_bayes_At0323.jpg", 
       width = 10, height = 3)




ggarrange(simu_rootn_bias_correct_density, simu_n_variance_correct_density, nrow = 1) -> longT_correct_density
annotate_figure(
  longT_correct_density # ,
  # top = text_grob("Density:correct", face = "bold", size = 14)
) ->  plot_longT_correct_density
plot_longT_correct_density

ggsave(plot_longT_correct_density, file="output/plots/plot_longT_correct_density0323.jpg", 
       width = 10, height = 3)







longT_correct_density_and_bayes <- ggarrange(
  ggarrange(
    text_grob(expression(psi[1]^"+" * "(" * hat(Q) * ")"), rot = 90, size = 12),
    plot_longT_correct_density,
    ncol = 2,
    widths = c(0.06, 1)
  ),
  ggarrange(
    text_grob(expression(psi[2]^"+" * "(" * hat(Q) * ")"), rot = 90, size = 12),
    plot_longT_correct_bayes,
    ncol = 2,
    widths = c(0.06, 1)
  ),
  nrow = 2,
  ncol = 1
)

ggsave(longT_correct_density_and_bayes, file="output/plots/plot_longT_correct_density_and_bayes0323.jpg", 
       width = 10, height = 3)
