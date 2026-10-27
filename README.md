# Inference for Causal Path-Specific Effects in Longitudinal Studies with Repeated Treatments and Multiple Mediators

🎯 This repository contains R analysis for paper "Inference for Causal Path-Specific Effects in Longitudinal Studies with Repeated Treatments and Multiple Mediators" ([arXiv](https://arxiv.org/abs/2610.34349)).

🎯 The code covers three simulation studies and a Framingham Heart Study (FHS) application examining the effects of sustained smoking on systolic blood pressure (SBP) through repeated body mass index (BMI) measurements. The simulations implement both density-based and treatment-regression one-step estimators. The FHS analysis uses the accompanying R package [**flexPaths**](https://github.com/xxou/flexPaths).

## Project overview

The simulation studies evaluate:

1.  **Correct nuisance-model specification:** finite-sample bias, variance, mean squared error, and coverage of the two one-step estimators.
2.  **Misspecified GLMs versus Super Learner:** performance under parametric model misspecification and flexible model estimation.
3.  **Robustness to model misspecification:** the multiple-robustness regimes considered in the manuscript, together with all-correct and all-misspecified benchmarks.

All three studies use two binary treatments and four ordered intermediate variables, with the temporal ordering `T1 → M1 → L1 → T2 → M2 → L2 → Y`. The simulation code uses `T1` and `T2` for the treatments denoted by $A_1$ and $A_2$ in the manuscript. The labels `density` and `bayes` refer to the density-based representation ($\ell=1$) and the Bayes-rule treatment-regression representation ($\ell=2$), respectively.

The FHS application uses smoking and BMI measurements from Exams **2, 4, 5, and 7**, with systolic blood pressure at Exam 7 as the outcome. It estimates the total smoking contrast, the contrast along selected BMI pathways, and the complementary pathway contrast.

## Folder description

``` text
longitudinal_PSE/
├── README.md
├── LICENSE.txt
├── longitudinal_PSE.Rproj
├── codes/
│   ├── functions/
│   ├── simulations/
│   │   ├── Theoretical_flexEst_0223/
│   │   │   └── DGP_coefficients/
│   │   └── misspecification_0222/
│   │       └── DGP_coefficients/
│   └── codes_FHS/
├── data/
│   └── FHS_Data.txt
└── output/
    ├── hpc_data/
    ├── local_parallel/
    │   ├── correct/
    │   └── miss/
    └── plots/
```

### `codes/functions/`

- [`est_funs_GLM.R`](codes/functions/est_funs_GLM.R) fits treatment models, conditional mediator models, treatment-history regressions, and sequential outcome regressions using GLMs. It constructs the nuisance quantities required by the simulation estimators.
- [`est_funs_SL.R`](codes/functions/est_funs_SL.R) constructs the corresponding nuisance quantities using Super Learner.
- [`one_step_estimator.R`](codes/functions/one_step_estimator.R) assembles the estimator components, computes one-step estimates and influence-function-based standard errors, and combines correctly and incorrectly specified nuisance estimates for the robustness simulations.
- [`SL.learner.R`](codes/functions/SL.learner.R) defines the custom `SL.dbarts` learner and prediction method, and the `SL.ranger.1000` wrapper with 1,000 trees.

These functions are written for the two-treatment simulation structure. The FHS application uses `flexPaths` to accommodate its four treatment occasions.

### `codes/simulations/Theoretical_flexEst_0223/`

This folder contains **Simulations 1 and 2**.

- [`1_longT_make_data.R`](codes/simulations/Theoretical_flexEst_0223/1_longT_make_data.R) defines the data generator and a function for Monte Carlo evaluation of counterfactual means.
- [`DGP_coefficients/par_0323.R`](codes/simulations/Theoretical_flexEst_0223/DGP_coefficients/par_0323.R) specifies the data-generating coefficients.
- [`3_simualtion_theoretical.R`](codes/simulations/Theoretical_flexEst_0223/3_simualtion_theoretical.R) runs Simulation 1 locally in parallel, summarizes the results, and produces the correct-specification figures. The filename retains its original spelling.
- [`2_complex_simulation.R`](codes/simulations/Theoretical_flexEst_0223/2_complex_simulation.R) contains the job-generation, single-task, aggregation, and plotting sections for Simulation 2. Its simulation stages are intended for execution on an HPC cluster; see the instructions below.

### `codes/simulations/misspecification_0222/`

This folder contains **Simulation 3**.

- [`1_makeData.R`](codes/simulations/misspecification_0222/1_makeData.R) defines the data generator and counterfactual simulation function for the robustness study.
- [`DGP_coefficients/par_0315.R`](codes/simulations/misspecification_0222/DGP_coefficients/par_0315.R) specifies its data-generating coefficients.
- [`2_simu_miss_falseFormula.R`](codes/simulations/misspecification_0222/2_simu_miss_falseFormula.R) specifies the nuisance-model configurations, runs the simulations locally in parallel, summarizes the estimates, and produces the robustness figures.

### `codes/codes_FHS/`

- [`0_clean_data2457.R`](codes/codes_FHS/0_clean_data2457.R) imports the FHS source files, derives variables, applies the missing-value updates and cohort selection rules, and saves the analysis datasets.
- [`variables.R`](codes/codes_FHS/variables.R) defines the baseline covariates, treatments, BMI mediators, time-varying covariate blocks, and outcome used in the analysis.
- [`1_descriptive.R`](codes/codes_FHS/1_descriptive.R) produces descriptive summaries and prints the baseline-characteristics table as LaTeX.
- [`2_PSE_task2457.R`](codes/codes_FHS/2_PSE_task2457.R) fits nuisance models with `flexPaths`, estimates the three counterfactual means, and computes their contrasts and confidence intervals. The active configuration uses Super Learner; saved numerical summaries for both Super Learner and GLM analyses appear in the final comments.

### `data/`

[`FHS_Data.txt`](data/FHS_Data.txt) provides the FHS data-access link. Individual-level FHS source data and the derived analysis datasets are **not included**. They must be obtained and generated before running the application.

### `output/`

- `hpc_data/` contains aggregated simulation results transferred from the HPC cluster. The current Simulation 2 plotting section reads [`longT.summarise.flex.SL.0323.Rdata`](output/hpc_data/longT.summarise.flex.SL.0323.Rdata).
- `local_parallel/correct/` contains saved correct-specification simulation results. The current Simulation 1 script reads and writes [`correct_specify0323_At.Rdata`](output/local_parallel/correct/correct_specify0323_At.Rdata).
- `local_parallel/miss/` contains saved robustness simulation results. The current Simulation 3 script writes [`GLM_misspecification0315.Rdata`](output/local_parallel/miss/GLM_misspecification0315.Rdata).
- [`plots/`](output/plots/) contains the simulation figures. The `0323` figures correspond to Simulations 1 and 2; the `0315_At` figures correspond to Simulation 3. `plot_misspecified.png` combines the robustness figures.

Other result files in these directories retain earlier simulation runs. Use the filenames listed above when following the current scripts.

## Citation and license

If you use this code, please cite [the paper](https://arxiv.org/abs/2610.34349):

``` bibtex
@article{ou2026causalpse,
  title   = {Inference for Causal Path-Specific Effects in Longitudinal Studies with Repeated Treatments and Multiple Mediators},
  author  = {Ou, Xiaxian and He, Xinwei and Nabi, Razieh},
  journal = {arXiv preprint arXiv:2610.34349},
  year    = {2026}
}
```

The code is distributed under the [MIT License](LICENSE.txt). Access to FHS data is governed separately by the study's data-access terms.
