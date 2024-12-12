# Performance and Stability Analysis of Winter Wheat Cultivars to Stagonospora nodorum Blotch Epidemics in Multi-environment Trials 

## Citation
The following source code accompanies our publication, currently under review in *Phytopathology* Journal:

```
@article{garnica202x,
  author = {Vinicius C. Garnica and Mohammad N. Shalizi and Peter S. Ojiambo},
  title = {Performance and stability of winter wheat cultivars to Stagonospora nodorum blotch epidemics in multi-environment trials},
  year = {202x},
  doi = {xxx},
  journal = {xxx}
}
```

## Introduction

This repository contains data and code for analyzing genotype-by-environment interaction (GEI) effects driving Stagonospora nodorum blotch (SNB) epidemics in winter wheat, based on small-plot research trials conducted in North Carolina. The study evaluates GEI effects using various variance-covariance models and key disease response metrics, including final foliar severity (**sev**), rate of incidence increase (**omega**), relative area under the disease progress stairs (**rAUDPS**), and time to 50% incidence (**T50**). Both genetic and non-genetic parameters are estimated using a third-order factor-analytic model, which provides insights into the performance and stability of winter wheat cultivars under varying environmental conditions.

This repository also outlines the complete analysis pipeline, detailing the necessary R scripts and data files for reproducing the study. The folder structure is organized as follows:

```
SNB_stability/
├── code/
│   ├── 0_stability_parameters.R
│   ├── 1_stability_analysis.R
│   ├── 2_factor_analytic_tools.R
│   ├── 3_plots.R
│   └── 4_NC_ecoregions_plot.R
├── data/
│   ├── raw_stability.Rdata
│   ├── res_fa_tools.Rdata
│   ├── res_stability.Rdata
│   ├── stability_final.Rdata
│   ├── stability_final.csv
│   └── map/
└── results/
│       ├── plots/
│       ├── CV.csv
│       ├── data_summary.csv
│       ├── exp_var.csv
│       ├── loadings.csv
│       ├── loadings.Rdata
│       ├── eblup_L1.csv
│       ├── op_rsmd.csv
│       ├── loadings.Rdata
└────── └── model_selection.csv

```

## Analysis pipeline

### Step 0: Stability Parameters

The `0_stability_parameters.R` script processes and analyzes data to assess cultivar stability through epidemiological metrics and disease progression models using the `raw_stability.RData` file. First, **sev** is created using the last disease assessment for each experimental unit. Second, **rAUDPS** is generated using the audps function from agricolae R package. Third, three non-linear population growth models (Monomolecular, Logistic, and Gompertz) are fitted to disease incidence data. Both **omega** and **T_50** are estimated from the best fitting model at each experimental unit. Results are saved on the file `stability_final.RData` for use in the next step.

### Step 1: Stability Analysis

The `1_stability_analysis.R` script processes the data in the `stability_final.RData` file, fitting various models and exporting the results to the `res_stability.RData` file. The output is a list containing the fits for the third-order factor analytic model. Within the script, you can adjust the model selection criteria as indicated in the code comments. 

### Step 2: Factor Analytic Tools

The `2_factor_analytic_tools.R` script applies the factor analytic selection tools (FAST) based on the previsouly selected third-order model and estimate genetic and non-genetic parameters. The fa.out function is a modified version of Smith et al. (2015) and Chaves et al. (2023). It exports the `res_fa_tools.RData` and `loadings.RData`. We use `loadings.RData` for **sev** in another study.

### Step 3: Plotting Results

Results are visualized using the `3_plots.R` script, which generates plots for different aspects of the analysis and associated csv files. `4_NC_ecoregions_plot.R` is for creating the ecoregion map of North Carolina. 

## References

* Smith, A.B., Ganesalingam, A., Kuchel, H. and Cullis, B.R., 2015. Factor analytic mixed models for the provision of grower information from national crop variety testing programs. Theoretical and applied genetics, 128, pp.55-72.

* Chaves, S.F., Evangelista, J.S., Trindade, R.S., Dias, L.A., Guimarães, P.E., Guimarães, L.J., Alves, R.S., Bhering, L.L. and Dias, K.O., 2023. Employing factor analytic tools for selecting high‐performance and stable tropical maize hybrids. Crop Science, 63(3), pp.1114-1125.



