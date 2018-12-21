# A nonparametric super-efficient estimator of the average treatment effect

**Authors:** David Benkeser, Wilson Cai, Mark van der Laan

-----

## Data analysis

This directory contains code and data needed to execute the analysis of the Cebuano Longitudinal Health Survey. The main script to execute is `analysis.R`, which sources in functions for the other `.R` scripts, and reads the `.RData` file. This code takes some time to execute owing to the need to cross-validate the super learner to obtain cross-validated standard error estimates. 