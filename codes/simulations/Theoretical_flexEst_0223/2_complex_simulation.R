# ************************************************************************************
# -----------------------            HPC: run                 -----------------------
# ************************************************************************************
# longT_complex0223
# use functions not package

# ....................................................
## --- (1) simulation: Rscript generate_job.R   ------
# ....................................................
unit_n = 1000

param_grid <- expand.grid(seed = 1:unit_n,
                          n = c(500, 1000, 2000, 4000, 8000),
                          model = c("GLM", "SL")
                          )

param_grid$joblist<-paste("Rscript task.R",
                          param_grid$seed,
                          param_grid$n,
                          param_grid$model)


n_list <- nrow(param_grid) / unit_n  # Write each chunk of simu_n rows to a separate text file

for (i in 1:n_list) {
  # Calculate the row indices for the ith chunk
  start_row <- ((i - 1) * unit_n) + 1
  end_row <- i * unit_n
  
  # Write the chunk to a text file
  write.table(param_grid[start_row:end_row, 'joblist'], 
              file = paste0("joblist_n",i,".txt"), 
              row.names = FALSE, col.names = FALSE, quote = F)
}



# ....................................................
## ------  (2) simulation: task.R              -------
# ....................................................

source("codes/simulations/Theoretical_flexEst_0223/1_longT_make_data.R")
source("codes/functions/est_funs_GLM.R")
source("codes/functions/est_funs_SL.R")
source("codes/functions/SL.learner.R")
source("codes/functions/one_step_estimator.R")
source("codes/simulations/Theoretical_flexEst_0223/DGP_coefficients/par_0323.R") 


pacman::p_load(tidyverse, dbarts, kernlab, SuperLearner, xgboost,earth, ranger)

args = commandArgs(trailingOnly=T)
seed = as.numeric(args[1])
n    = as.numeric(args[2])
model = as.character(args[3])


set.seed(seed) 
input_data = longT_make_data(n, par_T1, par_M1, par_L1, par_T2, par_M2, par_L2, par_Y)


# SL

if(model == "GLM"){
  
  # formula
  f.pro.T.list = list(pro.T1_formula = T1 ~ X1+X2+X3, 
                      pro.T2_formula = T2 ~ X1+X2+X3+T1+M1+L1)
  
  f.mediator.list = list(M1_formula = M1~X1+X2+X3+T1,
                         L1_formula = L1~X1+X2+X3+T1+M1,
                         M2_formula = M2~X1+X2+X3+T1+M1+L1+T2,
                         L2_formula = L2~X1+X2+X3+T1+M1+L1+T2+M2
  )
  
  f.bayes.list = list(T1_M1_formula = T1~X1+X2+X3+M1,
                      T1_L1_formula = T1~X1+X2+X3+M1+L1,
                      
                      T1_M2_formula = T1 ~ X1+X2+X3+M1+L1+T2+M2,
                      T1_L2_formula = T1 ~ X1+X2+X3+M1+L1+T2+M2+L2,
                      
                      T2_M2_formula = T2 ~ X1+X2+X3+M1+L1+M2+T1,
                      T2_L2_formula = T2 ~ X1+X2+X3+M1+L1+M2+L2+T1
  )
  
  seq.reg.list <- list(
    # E(Y| ....):
    u = predY ~ X1+X2+X3+M1+L1+M2+L2+T1+T2,
    # B1; consider E(L2|...) in addition:
    B1 = predY ~ X1+X2+X3+M1+L1+M2+T1+T2,
    # B2; consider E(M2|...) in addition:
    B2 = predY ~ X1+X2+X3+M1+L1+T1+T2,
    # B3; consider E(L1|...) in addition:
    B3 = predY ~ X1+X2+X3+M1+T1,
    # B4; consider E(M1|...) in addition
    B4 = predY ~ X1+X2+X3+T1
  )  
  
  
  funcs <- GLM_nuisance_funs(input_data, f.pro.T.list, f.mediator.list, f.bayes.list, seq.reg.list,
                             joint = "OR", At_correction = TRUE)
}


if(model == "SL"){
  
  SL_library = c('SL.xgboost','SL.ranger.1000', 'SL.dbarts', 'SL.glm', 'SL.glm.interaction') 
  
  funcs <- SL_nuisance_funs(data = input_data, X = c('X1', 'X2', 'X3'), SL_library = SL_library, joint = "OR")
}


component <- component.T1T2(input_data, funcs)
est.out = bind_rows(longT_one_step(component, type = "density")$out, 
                    longT_one_step(component, type = "bayes")$out)

est.out = est.out |> mutate(sample_size = n, model = model, seed = seed)


save(est.out, file = paste0("output/", model, ".output_", "n", n, "_", seed, ".Rdata" ))


# ...............................................................
## ------  (3) simulation: Rscript summarize.R           -------
# ...............................................................

param_grid <- expand.grid(seed = 1:1000,
                          n = c(500, 1000, 2000, 4000, 8000),
                          model = c("GLM", "SL")
)

param_grid$rerun_name <-paste("Rscript task.R",
                              param_grid$seed,
                              param_grid$n,
                              param_grid$model)

param_grid$file_name <-paste0("output/", 
                             param_grid$model, ".output_", 
                             "n", param_grid$n, "_", 
                             param_grid$seed, ".Rdata" )



false_index = c()
longT.summarise.flex.SL.0323 <- list()
for (i in 1:dim(param_grid)[1]) {
  # same as task.R code             
  file_name=param_grid$file_name[i]           
  
  # Attempt to load the data file      
  load_status <- tryCatch({        
    load(file_name)        
    TRUE       
  }, error = function(e) {
    FALSE     
  })       
  
  if(load_status==F){
    false_index = c(false_index,i)
  }
  
  if (load_status) {
    longT.summarise.flex.SL.0323[[i]] = est.out
  }
} 

save(longT.summarise.flex.SL.0323, file = "output/longT.summarise.flex.SL.0323.Rdata")

write.table(unique(param_grid[false_index, 'rerun_name']), 
            file = "rerun1.txt", 
            row.names = FALSE, col.names = FALSE, quote = F)



# **************************************************************************************
# ----------------------       Summary: results             ----------------------------
# **************************************************************************************
freshr::freshr()
pacman::p_load(tidyverse,gtsummary,gt, ggpubr)

load("output/hpc_data/longT.summarise.flex.SL.0323.Rdata")

simu.longT <- longT.summarise.flex.SL.0323 |>  map_df(~ .x)

# summary
simu.longT.summary <- simu.longT |> 
  mutate( true_value = 0.3935) |> 
  mutate(
    error = value - true_value,
    square_error = (value - true_value)^2,
    CI.upper = value + 1.96*SE,
    CI.lower = value - 1.96*SE,
    cover = if_else( (true_value <= CI.upper) & (true_value >= CI.lower), 1, 0),
    length = CI.upper - CI.lower)

simu.longT.evaluate <- simu.longT.summary |>   
  group_by(type, sample_size, model) |>  
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
    
    true_var = 13.9
  )



### (1) display table 
gt.simu.longT.evaluate <- simu.longT.evaluate |> 
  select(type, model, sample_size, bias, sd_est, mse, cover_rate, avg_length) |> 
  pivot_wider(
    names_from = model,
    values_from = c(bias, sd_est, mse, cover_rate, avg_length),
    names_glue = "{.value}_{model}" 
  ) |> 
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
model_labels <- c(
  "SL" = "Super Learner",
  "glm" = "GLM"
)

# root n bias
simu.longT.evaluate.bayes <- simu.longT.evaluate |> filter(type == "bayes")
simu.longT.evaluate.density <- simu.longT.evaluate |> filter(type == "density")

ggplot(simu.longT.evaluate.bayes, aes(x = sample_size, y = scaled_bias)) +
  geom_point(color = "red", size = 2) +
  geom_line() +
  geom_errorbar(aes(ymin = CI.low.scaled_bias, ymax = CI.up.scaled_bias), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "blue") +
  labs(x = "Sample Size", y = expression(sqrt(n)-Bias)) +
  theme_bw() +
  facet_grid(rows = vars(model),scales = "fixed") +
  scale_x_continuous(
    breaks = unique(simu.longT.evaluate.bayes$sample_size),  # Set breaks at sample_size values
    labels = unique(simu.longT.evaluate.bayes$sample_size)  # Set labels as sample_size
  )->simu_rootn_bias_bayes


ggplot(simu.longT.evaluate.density, aes(x = sample_size, y = scaled_bias)) +
  geom_point(color = "red", size = 2) +
  geom_line() +
  geom_errorbar(aes(ymin = CI.low.scaled_bias, ymax = CI.up.scaled_bias), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "blue") +
  labs(x = "Sample Size", y = expression(sqrt(n)-Bias)) +
  theme_bw() +
  facet_grid(rows = vars(model),scales = "fixed") +
  scale_x_continuous(
    breaks = unique(simu.longT.evaluate.density$sample_size),  # Set breaks at sample_size values
    labels = unique(simu.longT.evaluate.density$sample_size)  # Set labels as sample_size
  )->simu_rootn_bias_density



# root n variance
ggplot(simu.longT.evaluate.bayes, aes(x = sample_size, y = scaled_var_est)) +
  geom_point(color = "red", size = 2) +
  geom_line() +
  geom_errorbar(aes(ymin = CI.low.scaled_var_est, ymax = CI.up.scaled_var_est), width = 0.2) +
  labs(x = "Sample Size", y = expression(n-Variance)) +
  theme_bw() +
  geom_hline(yintercept = simu.longT.evaluate.bayes$true_var, linetype = "dashed", color = "blue") +
  facet_grid(rows = vars(model), scales = "fixed") +
  scale_x_continuous(
    breaks = unique(simu.longT.evaluate.bayes$sample_size),  # Set breaks at sample_size values
    labels = unique(simu.longT.evaluate.bayes$sample_size)  # Set labels as sample_size
  )->simu_n_variance_bayes


ggplot(simu.longT.evaluate.density, aes(x = sample_size, y = scaled_var_est)) +
  geom_point(color = "red", size = 2) +
  geom_line() +
  geom_errorbar(aes(ymin = CI.low.scaled_var_est, ymax = CI.up.scaled_var_est), width = 0.2) +
  labs(x = "Sample Size", y = expression(n-Variance)) +
  theme_bw() +
  geom_hline(yintercept = simu.longT.evaluate.bayes$true_var, linetype = "dashed", color = "blue") +
  facet_grid(rows = vars(model), scales = "fixed") +
  scale_x_continuous(
    breaks = unique(simu.longT.evaluate.density$sample_size),  # Set breaks at sample_size values
    labels = unique(simu.longT.evaluate.density$sample_size)  # Set labels as sample_size
  )->simu_n_variance_density


ggarrange(simu_rootn_bias_bayes, simu_n_variance_bayes, nrow = 1) -> longT_flex.SL_bayes
annotate_figure(
  longT_flex.SL_bayes #,
#  top = text_grob("Bayes:SL.xgboost, SL.ranger.1000, SL.glm, SL.dbarts, SL.glm.interaction\n", face = "bold", size = 14)
) ->  plot_longT_flex.SL_bayes
plot_longT_flex.SL_bayes

ggsave(plot_longT_flex.SL_bayes, file="output/plots/plot_longT_flex.SL_bayes0323.jpg", 
       width = 10, height = 4)




ggarrange(simu_rootn_bias_density, simu_n_variance_density, nrow = 1) -> longT_flex.SL_density
annotate_figure(
  longT_flex.SL_density #,
 # top = text_grob("Density:SL.xgboost, SL.ranger.1000, SL.glm, SL.dbarts, SL.glm.interaction\n", face = "bold", size = 14)
) ->  plot_longT_flex.SL_density
plot_longT_flex.SL_density


ggsave(plot_longT_flex.SL_density, file="output/plots/plot_longT_flex.SL_density0323.jpg", 
       width = 10, height = 4)










