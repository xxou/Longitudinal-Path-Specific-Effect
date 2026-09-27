
# ************************************************************************************
# -----------------------            Simulation                -----------------------
## Feb - 22 - 2025
# Simplest data for misspecification
# GLM : locally
# False formula
# ************************************************************************************
# longT_complex0223
# use functions not package

freshr::freshr()

# locally 
pacman::p_load(tidyverse, matrixStats)


source("codes/simulations/misspecification_0222/1_makeData.R")
source("codes/functions/est_funs_GLM.R")
source("codes/functions/one_step_estimator.R")


# **************************************
### one simulation -----------------

run_simulation <- function(seed, n, type = "both"){
  
  set.seed(seed)
  cat("seed:", seed, "; n:",  n, "\n")
  
  
  
  # *********  Generate data  ***************
  source("codes/simulations/misspecification_0222/DGP_coefficients/par_0315.R")
  
  data = longT_make_simplest_data(n = n, par_T1, par_M1, par_L1, par_T2, par_M2, par_L2, par_Y)
  
  
  
  
  # ********* correct ***************
  
  correct_funcs <- GLM_nuisance_funs(data, 
                                     f.pro.T.list = list(pro.T1_formula = T1 ~ X1+X2+X3+X4, 
                                                         pro.T2_formula = T2 ~ X1+X2+X3+X4+T1+M1+L1), 
                                     
                                     f.mediator.list = list(M1_formula = M1~X1+X2+X3+X4+T1,
                                                            L1_formula = L1~X1+X2+X3+X4+T1+M1,
                                                            M2_formula = M2~X1+X2+X3+X4+T1+M1+L1+T2,
                                                            L2_formula = L2~X1+X2+X3+X4+T1+M1+L1+T2+M2), 
                                     
                                     f.bayes.list = list(T1_M1_formula = T1~X1+X2+X3+X4+M1,
                                                         T1_L1_formula = T1~X1+X2+X3+X4+M1+L1,
                                                         
                                                         T1_M2_formula = T1~X1+X2+X3+X4+M1+L1+T2+M2,
                                                         T1_L2_formula = T1~X1+X2+X3+X4+M1+L1+T2+M2+L2,
                                                         
                                                         T2_M2_formula = T2~X1+X2+X3+X4+T1+M1+L1+M2,
                                                         T2_L2_formula = T2~X1+X2+X3+X4+T1+M1+L1+M2+L2), 
                                     
                                     seq.reg.list = list(
                                       # E(Y| ....):
                                       u = predY ~ X1+X2+X3+X4+M1+L1+M2+L2+T1+T2,
                                       # B1; consider E(L2|...) in addition:
                                       B1 = predY ~ X1+X2+X3+X4+M1+L1+M2+T1+T2,
                                       # B2; consider E(M2|...) in addition:
                                       B2 = predY ~ X1+X2+X3+X4+M1+L1+T1+T2,
                                       # B3; consider E(L1|...) in addition:
                                       B3 = predY ~ X1+X2+X3+X4+M1+T1,
                                       # B4; consider E(M1|...) in addition
                                       B4 = predY ~ X1+X2+X3+X4+T1
                                     ),
                                     
                                     joint = "OR", At_correction = TRUE
  )
  
  
  correct_component <- component.T1T2(data,  correct_funcs)
  correct_est <- rbind(longT_one_step( correct_component, type = "bayes")$out, 
                       longT_one_step( correct_component, type = "density")$out)
  correct_est$condition <- "correct"
  
  
  
  # ********* false ***************
  
  false_funcs <- GLM_nuisance_funs(data, 
                                   f.pro.T.list = list(pro.T1_formula = T1 ~ I(X1^2)+I(exp(X2))+I(X1^0.5/X3)+I(cos(10*X4)), 
                                                       pro.T2_formula = T2 ~ I(X1^2)+I(exp(X2))+I(X1^0.5/X3)+I(cos(10*X4))+T1+M1+L1), 
                                   
                                   f.mediator.list = list(M1_formula = M1~I(X1^2)+I(exp(X2))+I(X1^0.5/X3)+I(cos(10*X4))+T1,
                                                          L1_formula = L1~I(X1^2)+I(exp(X2))+I(X1^0.5/X3)+I(cos(10*X4))+T1+M1,
                                                          M2_formula = M2~I(X1^2)+I(exp(X2))+I(X1^0.5/X3)+I(cos(10*X4))+T1+M1+L1+T2,
                                                          L2_formula = L2~I(X1^2)+I(exp(X2))+I(X1^0.5/X3)+I(cos(10*X4))+T1+M1+L1+T2+M2), 
                                   
                                   f.bayes.list = list(T1_M1_formula = T1~I(X1^2)+I(exp(X2))+I(X1^0.5/X3)+I(cos(10*X4))+M1,
                                                       T1_L1_formula = T1~I(X1^2)+I(exp(X2))+I(X1^0.5/X3)+I(cos(10*X4))+M1+L1,
                                                       
                                                       T1_M2_formula = T1~I(X1^2)+I(exp(X2))+I(X1^0.5/X3)+I(cos(10*X4))+M1+L1+T2+M2,
                                                       T1_L2_formula = T1~I(X1^2)+I(exp(X2))+I(X1^0.5/X3)+I(cos(10*X4))+M1+L1+T2+M2+L2,
                                                       
                                                       T2_M2_formula = T2~I(X1^2)+I(exp(X2))+I(X1^0.5/X3)+I(cos(10*X4))+T1+M1+L1+M2,
                                                       T2_L2_formula = T2~I(X1^2)+I(exp(X2))+I(X1^0.5/X3)+I(cos(10*X4))+T1+M1+L1+M2+L2), 
                                   
                                   seq.reg.list = list(
                                     # E(Y| ....):
                                     u = predY ~ I(X1^2)+I(exp(X2))+I(X1^0.5/X3)+I(cos(10*X4))+T1+M1+L1+M2+L2+T2,
                                     # B1; consider E(L2|...) in addition:
                                     B1 = predY ~ I(X1^2)+I(exp(X2))+I(X1^0.5/X3)+I(cos(10*X4))+T1+M1+L1+M2+T2,
                                     # B2; consider E(M2|...) in addition:
                                     B2 = predY ~ I(X1^2)+I(exp(X2))+I(X1^0.5/X3)+I(cos(10*X4))+T1+M1+L1+T2,
                                     # B3; consider E(L1|...) in addition:
                                     B3 = predY ~ I(X1^2)+I(exp(X2))+I(X1^0.5/X3)+I(cos(10*X4))+T1+M1,
                                     # B4; consider E(M1|...) in addition
                                     B4 = predY ~ I(X1^2)+I(exp(X2))+I(X1^0.5/X3)+I(cos(10*X4))+T1
                                   ),
                                   
                                   joint = "OR", At_correction = TRUE
  )
  
  false_funcs <- purrr::map_if(
    false_funcs, names(false_funcs) %in% c("prob1.T1", "prob0.T1", "prob1.T2", "prob0.T2"),  ~ pmax(pmin(.x, 0.99), 0.01)  
  )
  
  
  false_component <- component.T1T2(data,  false_funcs)
  false_est <-rbind(longT_one_step( false_component, type = "density")$out,
                    longT_one_step( false_component, type = "bayes")$out) 
  false_est$condition <- "false"
  
  
  d.con1_est = d.con2_est = d.con3_est = b.con_est = NULL
  
  # ************************************************
  ##  Density ratio  -----------------
  # ************************************************
  
  if( (type == "both") | (type == "density") ){
  # ********* condition1 ***************
  # f(T1|..), f(T2|..), f(L1|..), f(L2|..)
  d.con1_funcs <- GLM_nuisance_funs(data, 
                                  f.pro.T.list = list(pro.T1_formula = T1 ~ X1+X2+X3+X4, 
                                                      pro.T2_formula = T2 ~ X1+X2+X3+X4+T1+M1+L1), 
                                  
                                  f.mediator.list = list(M1_formula = M1~I(X1^2)+I(exp(X2))+I(X1^0.5/X3)+I(cos(10*X4))+T1,
                                                         L1_formula = L1~X1+X2+X3+X4+T1+M1,
                                                         M2_formula = M2~I(X1^2)+I(exp(X2))+I(X1^0.5/X3)+I(cos(10*X4))+T1+M1+L1+T2,
                                                         L2_formula = L2~X1+X2+X3+X4+T1+M1+L1+T2+M2), 
                                  f.bayes.list = NULL, 
                                  
                                  seq.reg.list = list(
                                    # E(Y| ....):
                                    u = predY ~ I(X1^2)+I(exp(X2))+I(X1^0.5/X3)+I(cos(10*X4))+M1+L1+M2+L2+T1+T2,
                                    
                                    # B1; consider E(L2|...) in addition:
                                    B1 = predY ~ I(X1^2)+I(exp(X2))+I(X1^0.5/X3)+I(cos(10*X4))+X1+X2+X3+X4+M1+L1+M2+T1+T2,
                                    # B1 = predY ~ X1+X2+X3+X4+M1+L1+M2+T1+T2,
                                    
                                    # B2; consider E(M2|...) in addition:
                                    B2 = predY ~ I(X1^2)+I(exp(X2))+I(X1^0.5/X3)+I(cos(10*X4))+X1+X2+X3+X4+M1+L1+T1+T2,
                                    
                                    # B3; consider E(L1|...) in addition:
                                    B3 = predY ~ I(X1^2)+I(exp(X2))+I(X1^0.5/X3)+I(cos(10*X4))+X1+X2+X3+X4+M1+T1,
                                    # B3 = predY ~ X1+X2+X3+X4+M1+T1,
                                    
                                    # B4; consider E(M1|...) in addition
                                    B4 = predY ~ I(X1^2)+I(exp(X2))+I(X1^0.5/X3)+I(cos(10*X4))+X1+X2+X3+X4+T1) 
                                  )
  
  
  d.con1_component <- component.T1T2(data, d.con1_funcs)
  d.con1_est <- longT_one_step(d.con1_component, type = "density")$out
  d.con1_est$condition <- "condition 1"
  
  
  # ********* condition2 ***************
  # outcome regression E(Y|..), f(M1|..), f(M2|..), f(L1|..), f(L2|..)
  
  d.con2_funcs <- GLM_nuisance_funs(data, 
                                  f.pro.T.list = list(pro.T1_formula = T1 ~ I(X1^2)+I(exp(X2))+I(X1^0.5/X3)+I(cos(10*X4)), 
                                                      pro.T2_formula = T2 ~ I(X1^2)+I(exp(X2))+I(X1^0.5/X3)+I(cos(10*X4))+T1+M1+L1), 
                                  
                                  f.mediator.list = list(M1_formula = M1~X1+X2+X3+X4+T1,
                                                         L1_formula = L1~X1+X2+X3+X4+T1+M1,
                                                         M2_formula = M2~X1+X2+X3+X4+T1+M1+L1+T2,
                                                         L2_formula = L2~X1+X2+X3+X4+T1+M1+L1+T2+M2), 
                                  
                                  f.bayes.list = NULL, 
                                  
                                  seq.reg.list = list(
                                    # E(Y| ....):
                                    u = predY ~ X1+X2+X3+X4+M1+L1+M2+L2+T1+T2,
                                    # B1; consider E(L2|...) in addition:
                                    B1 = predY ~ X1+X2+X3+X4+M1+L1+M2+T1+T2,
                                    # B2; consider E(M2|...) in addition:
                                    B2 = predY ~ X1+X2+X3+X4+M1+L1+T1+T2,
                                    # B3; consider E(L1|...) in addition:
                                    B3 = predY ~ X1+X2+X3+X4+M1+T1,
                                    # B4; consider E(M1|...) in addition
                                    B4 = predY ~ X1+X2+X3+X4+T1) 
  )
  
  d.con2_funcs <- purrr::map_if(
    d.con2_funcs, names(d.con2_funcs) %in% c("prob1.T1", "prob0.T1", "prob1.T2", "prob0.T2"),  ~ pmax(pmin(.x, 0.99), 0.01)  
  )
  

  d.con2_component <- component.T1T2(data, d.con2_funcs)
  d.con2_est <- longT_one_step(d.con2_component, type = "density")$out
  d.con2_est$condition <- "condition 2"
  
  
  
  # ********* condition3 ***************
  # f(T1|..), f(T2|..), E(Y|..), f(M1|..), f(M2|..)
  
  d.con3_funcs <- GLM_nuisance_funs(data, 
                                  f.pro.T.list = list(pro.T1_formula = T1 ~ X1+X2+X3+X4, 
                                                      pro.T2_formula = T2 ~ X1+X2+X3+X4+T1+M1+L1), 
                                  
                                  f.mediator.list = list(M1_formula = M1~X1+X2+X3+X4+T1,
                                                         L1_formula = L1~I(X1^2)+I(exp(X2))+I(X1^0.5/X3)+I(cos(10*X4))+T1+M1,
                                                         M2_formula = M2~X1+X2+X3+X4+T1+M1+L1+T2,
                                                         L2_formula = L2~I(X1^2)+I(exp(X2))+I(X1^0.5/X3)+I(cos(10*X4))+T1+M1+L1+T2+M2), 
                                  
                                  f.bayes.list = NULL, 
                                  
                                  seq.reg.list = list(
                                    # E(Y| ....):
                                    u = predY ~ X1+X2+X3+X4+M1+L1+M2+L2+T1+T2,
                                    # B1; consider E(L2|...) in addition:
                                    B1 = predY ~ X1+X2+X3+X4+I(X1^2)+I(exp(X2))+I(X1^0.5/X3)+I(cos(10*X4))+M1+L1+M2+T1+T2,
                                    # B2; consider E(M2|...) in addition:
                                    B2 = predY ~ X1+X2+X3+X4+I(X1^2)+I(exp(X2))+I(X1^0.5/X3)+I(cos(10*X4))+M1+L1+T1+T2,
                                    # B3; consider E(L1|...) in addition:
                                    B3 = predY ~ X1+X2+X3+X4+I(X1^2)+I(exp(X2))+I(X1^0.5/X3)+I(cos(10*X4))+M1+T1,
                                    # B4; consider E(M1|...) in addition
                                    B4 = predY ~ X1+X2+X3+X4+I(X1^2)+I(exp(X2))+I(X1^0.5/X3)+I(cos(10*X4))+T1) 
                                    
  )
  
  
  d.con3_component <- component.T1T2(data, d.con3_funcs)
  d.con3_est <- longT_one_step(d.con3_component, type = "density")$out
  d.con3_est$condition <- "condition 3"
  
  }
  
  
  # ************************************************
  ##  Bayes ratio  -----------------
  # ************************************************
  
  if( (type == "both") | (type == "bayes") ){
  # use false data strategy
  b.con_funcs <- condtion.funcs(correct_funcs, false_funcs, type = "bayes")
  b.con_est <- map(names(b.con_funcs), ~{
    
    component <- component.T1T2(data, b.con_funcs[[.x]])
    result <- longT_one_step(component, type = "bayes")$out
    result$condition <- paste0("condition ", str_extract(.x, "\\d+"))
    
    return(result)
  })
  
  }
  # *********  est.out  ***************
  est.out = bind_rows(d.con1_est, d.con2_est, d.con3_est, b.con_est, correct_est, false_est)
  
  est.out = est.out |> mutate(sample_size = n,  seed = seed)
  
  return(est.out)
}

run_simulation(seed = 153, n=2000)


# **************************************
### parallel simulation -----------------

library(furrr)

plan(multisession)
# setting  ********************
unit_n = 1000

param_grid <- expand.grid(seed = 1:unit_n,
                          n = c(500, 1000, 2000, 4000, 8000)  # 500, 1000, 2000, 
)

nbreaks = 20
results = vector('list', nbreaks)
param_grid_splits <- split(param_grid, cut(seq_len(nrow(param_grid)), breaks = nbreaks, labels = FALSE))

for (i in seq_len(nbreaks)) {
  
  cat("i:", i, " ------------------------ ", "\n")
  
  results[[i]] <- future_map2(param_grid_splits[[i]]$seed, param_grid_splits[[i]]$n,
                              run_simulation, .options = furrr_options(seed = NULL))
}

# Combine results into a dataframe
est.out <- bind_rows(results)
save(est.out, file =  "output/local_parallel/miss/GLM_misspecification0315.Rdata")

# GLM_misspecification0310_absC.Rdata:  B2 = predY ~ X1+X2+X3+X4+M1+L1+T1+T2
# GLM_misspecification0310_revC.Rdata:  B2 = predY ~ I(X1^2)+I(exp(X2))+I(X1^0.5/X3)+I(cos(10*X4))+X1+X2+X3+X4+M1+L1+T1+T2



# **************************************************************************************
# ----------------------       Summary: results             ----------------------------
# **************************************************************************************
# freshr::freshr()
pacman::p_load(tidyverse,gtsummary,gt, forcats,ggpubr,patchwork,grid,gridExtra)

quantile(est.out$value, c(0.00005,0.99995))
est.out |> arrange(desc(value)) |> head(12) 
est.out |> arrange((value)) |> mutate(value = round(value,3))|> head(20)

# est.out <- est.out |> group_by(type,condition, sample_size) |>  filter( (value>= quantile(value, 0.01)) &
#                                                                         (value <= quantile(value, 0.99))
#                                                                         )

est.out <- est.out |>   filter( (value>= 1) &
                                  (value <=10)
)

dim(est.out)

## summary output --------------
# load("output/local_parallel/miss/GLM_misspecification0310_absC.Rdata")



mis.simu.summary <- est.out |> 
  mutate(
    true_est = 2.201,
    true_var = 13.053
  ) |> 
  mutate(
    error = value - true_est,
    square_error = (value - true_est)^2,
    CI.upper = value+1.96*SE,
    CI.lower = value-1.96*SE,
    cover = if_else( (true_est <= CI.upper) & (true_est>=CI.lower),1,0),
    length = CI.upper-CI.lower)

head(mis.simu.summary)

mis.simu.evaluate <- mis.simu.summary |>   
  group_by(type,condition, sample_size) |>  
  mutate(dev = value-mean(value),
         square_dev = (value-mean(value))^2) |> 
  summarise(
    bias=mean(error),
    se_bias = sqrt(mean(square_dev)/n()),
    mse = mean(square_error),
    
    sd_est = sd(value),
    var_est = mean(square_dev),
    se_var_est = sqrt(var(square_dev)/n()),
    
    cover_rate = mean(cover,na.rm=T),
    avg_length =  mean(length,na.rm=T) ,
    
    true_var = mean(true_var)
    
  ) |>  ungroup() |>  
  mutate(
    scaled_bias = sqrt(sample_size)*bias,
    se_scaled_bias = sqrt(sample_size)*se_bias,
    CI.up.scaled_bias = scaled_bias+1.96*se_scaled_bias,
    CI.low.scaled_bias = scaled_bias-1.96*se_scaled_bias,
    
    scaled_var_est = sample_size*var_est,
    se_scaled_var_est = sample_size*se_var_est,
    CI.up.scaled_var_est = scaled_var_est+1.96*se_scaled_var_est,
    CI.low.scaled_var_est = scaled_var_est-1.96*se_scaled_var_est
  )




### (3) plot summary -----------

plot_list <- list()
Type.combination <- unique(mis.simu.evaluate$type)

for (i in seq_along(Type.combination)) {
  Type <- Type.combination[i]
  
  # Filter the data for the current Path
  data_type <- mis.simu.evaluate %>% filter(type == Type)
  data_type$condition <- fct_relevel(factor(data_type$condition), "correct", "false", after = Inf)
  
  # Plot for scaled_bias
  plot_bias <- ggplot(data_type, aes(x = sample_size, y = scaled_bias)) +
    geom_point(color = "red", size = 2) +
    geom_line() +
    geom_errorbar(aes(ymin = CI.low.scaled_bias, ymax = CI.up.scaled_bias), width = 0.2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "blue") +
    labs(x = "Sample Size", y = expression(sqrt(n) - Bias)) +
    theme_bw(base_size = 16) +  # Increase the base font size for the whole theme
    facet_grid(rows = vars(condition), scales = "fixed") +
    scale_x_continuous(
      breaks = unique(data_type$sample_size),  # Set breaks at sample_size values
      labels = unique(data_type$sample_size)  # Set labels as sample_size
    ) +
    theme(
      axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 12),
      axis.title = element_text(size = 14, face = "bold"),  # Increase axis label size
      strip.text = element_text(size = 14),  # Increase facet label size
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5)  # Adjust title size
    )
  
  # Plot for scaled_var_est
  plot_var <- ggplot(data_type, aes(x = sample_size, y = scaled_var_est)) +
    geom_point(color = "red", size = 2) +
    geom_line() +
    geom_errorbar(aes(ymin = CI.low.scaled_var_est, ymax = CI.up.scaled_var_est), width = 0.2) +
    geom_hline(yintercept = data_type$true_var[1], linetype = "dashed", color = "blue") +
    labs(x = "Sample Size", y = expression(n - Variance)) +
    theme_bw(base_size = 16) +  # Increase the base font size for the whole theme
    facet_grid(rows = vars(condition), scales = "fixed") +
    scale_x_continuous(
      breaks = unique(data_type$sample_size),  # Set breaks at sample_size values
      labels = unique(data_type$sample_size)  # Set labels as sample_size
    ) +
    theme(
      axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 12),
      axis.title = element_text(size = 14, face = "bold"),  # Increase axis label size
      strip.text = element_text(size = 14),  # Increase facet label size
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5)  # Adjust title size
    )
  
  # Combine the two plots using gridExtra
  combined_plot <- grid.arrange(
    arrangeGrob(
      plot_bias,
      plot_var,
      ncol = 2
      # ,
      # left = textGrob(
      #   latex_labels[i],
      #   rot = 90,  # Rotate to make it a vertical axis label
      #   gp = gpar(fontsize = 24, fontface = "bold")  # Larger font for y-axis label
      # )
    )
  )
  
  # Save the combined plot in the list
  plot_list[[Type]] <- combined_plot
  
  # Optional: Save the combined plot as a file
 ggsave(filename = paste0("output/plots/plot_misspecification_", Type, "0315_At.png"), plot = combined_plot, width = 15, height = 12)
}


whole_plot <- ggarrange(plotlist = plot_list[1:5], ncol = 1, nrow = 5)
ggsave(filename = paste0("output/plots/plot_misspecified.png"), plot = whole_plot, width = 14, height = 28)


