# Performance and Stability Analysis of Winter Wheat Cultivars to Stagonospora nodorum Blotch Epidemics in Multi-environment Trials 

## Citation
This is the source code for our publication::

```
@article{garnica202x,
  author = {Vinicius C. Garnica and Mohammad N. Shalizi and Peter S. Ojiambo},
  title = {Performance and Stability Analysis of Winter Wheat Cultivars to Stagonospora nodorum Blotch Epidemics in Multi-environment Trials},
  year = {202x},
  doi = {xxx},
  journal = {xxx}
}
```

## Introduction

This repository contains data and code for analyzing genetic-by-environment interaction (GEI) effects driving Stagonospora nodorum blotch (SNB) epidemics in wheat, based on small-plot research trials in North Carolina. Our study investigates GEI effects using various variance-covariance models and disease response metrics, including final foliar severity, rate of incidence increase, rAUDPS, and T50. We report both genetic and non-genetic parameters derived from a factor analytic model.

The pipeline used for the analysis, emphasizing the required Rdata files and R scripts for each step, is described below.

The folder structure and files are as follows:

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
│   └── map/
└── results/
│       ├── plots/
│       ├── CV.csv
│       ├── data_summary.csv
│       ├── exp_var.csv
│       ├── loadings.csv
│       ├── eblup_L1.csv
│       ├── op_rsmd
│       ├── loadings.Rdata
└────── └── model_selection.csv

 
```

## Analysis pipeline

### Step 0: Stability Parameters

The `0_stability_parameters.R` script creates responses used in the study using the `stability_final.Rdata` file. This script exports a Rdata file named `weather_data.Rdata` for use in the next step.

### Step 1: Stability Analysis

The `1_stability_analysis.R` script fits various models to the data and exports the `res_stability.RData` file, which is a list of multiple parameters and models. In the `1_stability_analysis.R` script, you may be able to change the model selection criteria, as noted in the code. 

### Step 2: Factor Analytic Tools

The `2_factor_analytic_tools.R` script applies the factor analytic model you selected to estimate genetic and non-genetic parameters. The fa.out function is a modified version of Smith et al. (2015) and Chaves et al. (2023). It exports the `res_fa_tools.RData` and `loadings.RData`. We use these loadings of SNB foliar severity in another study.

### Step 3: Plotting Results

Results are visualized using the `3_plots.R` script, which generates plots for different aspects of the analysis and associated csv files. `4_NC_ecoregions_plot.R` is for creating the ecoregion map of North Carolina. 

## References

* Smith, A.B., Ganesalingam, A., Kuchel, H., et al. (2015). Factor analytic mixed models for the provision of grower information from national crop variety testing programs. Theoretical and Applied Genetics, 128, 55–72.

* Chaves, S.F.S., Evangelista, J.S.P.C., Trindade, R.S., Dias, L.A.S., Guimarães, P.E., Guimarães, L.J.M., Alves, R.S., Bhering, L.L., & Dias, K.O.G. (2023). Employing factor analytic tools for selecting high-performance and stable tropical maize hybrids. Crop Science, 63, 1114–1125.



