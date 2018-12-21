# A nonparametric super-efficient estimator of the average treatment effect

**Authors:** David Benkeser, Wilson Cai, Mark van der Laan

-----

## Description 

This repository contains code to reproduce the simulation and data analysis in the paper 
"A nonparametric super-efficient estimator of the average treatment effect". All code was implemented in the `R` programming language. 

-----

## Directory

The `sim1` folder contains code that was used in the first simulation that was based on parametric models. The `sim2` folder contains code that was used in the second simulation study based on the famous Kang and Schafer simulation design. Both simulations rely on an `R` package developed by our colleagues Jeremy Coyle and Nima Hejazi, which implements the HAL MLE estimator. This package is in development and is not yet available on CRAN. A developmental version can be downloaded from GitHub as follows.

```r
# need devtools package
devtools::install_github("jeremyrcoyle/hal9001")
```

The `data_analysis` folder contains code that was used in the analysis of the Cebuano Longitudinal Health and Nutrition Survey (CLHNS). The raw data can be downloaded from a the [Cebu Longitudinal Health and Nutrition Survey Dataverse](https://dataverse.unc.edu/dataverse/cebu). We include a cleaned version of the data in `.RData` format. 

Each directory has its own README document to summarize the contents. 