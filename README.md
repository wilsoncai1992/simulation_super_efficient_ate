# Simulation of Benkeser D, Cai W, van der Laan MJ (2019+). *A nonparametric super-efficient estimator of the average treatment effect.*


## Pre-requisites

```R
devtools::install_github("tlverse/hal9001")
install.packages(c("glmnet", "dplyr", "drtmle", "tidyverse", "Rmpi", "doMPI", "SuperLearner", "ggpubr"))
```


## Instructions

To reproduce the simulation section

```R
# run simulation 1
cd ./code_iv/
R CMD BATCH test-oat_iv.R
mkdir ./output/
mv ./df_mc_result.rda ./output/df_mc_result.rda
# create plots
R CMD BATCH ./plot_results_paper.R
# the plots will be saved in `./code_iv/output/` after the script

# run simulation 2
cd ./code_ks/
R CMD BATCH test-oat_ks.R
mkdir ./output/
mv ./df_mc_result.rda ./output/df_mc_result.rda
# create plots
R CMD BATCH ./plot_results_paper.R
# the plots will be saved in `./code_ks/output/` after the script
```

To reproduce the real data analysis section

```R
# run real data analysis
cd ./code_real_data/catnap/
R CMD BATCH run_npbootstrap.R
R CMD BATCH run_permutey.R
mkdir ./result_npboot
mkdir ./result_permutey
mv npbootstrap.rda ./result_npboot/npbootstrap.rda
mv permute_y.rda ./result_permutey/permute_y.rda
# create plots
R CMD BATCH ./plot.R
# the plots will be saved in `./code_real_data/result_npboot/` and `./code_real_data/result_permutey/` after the script
```

