# BayesCOOP

The repository contains the analysis codes from the Bayesian coopertive
project.

## Dependencies

`BayesCOOP` requires the following `R` package: `devtools` (for
installation only). Please install it before installing `BayesCOOP`,
which can be done as follows (execute from within a fresh R session):

    install.packages("devtools")
    library(devtools)

## Installation

Once the dependencies are installed, `BayesCOOP` can be loaded using the
following command:

    devtools::install_github("himelmallick/BayesCOOP")
    library(BayesCOOP)

## Load libraries

    library(BhGLM)
    library(tidyverse)
    library(caret)

## Sourcing required functions

    source("~/bayesCOOP/BC/bayesCoop.R")

## Loading the StelzerDOS real dataset

    data_train = get(load(url("https://raw.githubusercontent.com/himelmallick/IntegratedLearner/master/data/StelzerDOS.RData"))); 
    rm(pcl)

    data_test = get(load(url("https://raw.githubusercontent.com/himelmallick/IntegratedLearner/master/data/StelzerDOS_valid.RData"))); 
    rm(pcl)

## Pre-processing the longitudinal data by considering only baseline observations

### Remove metabolomics from the train set to match with validation

    data_train$feature_metadata = data_train$feature_metadata %>% dplyr::filter(featureType!='Metabolomics')
    data_train$feature_table = data_train$feature_table[rownames(data_train$feature_metadata),]

### Consider only baseline observations for the train set

    positions = grep("A", colnames(data_train$feature_table), ignore.case = TRUE)
    data_train$feature_table = data_train$feature_table[, positions]
    data_train$sample_metadata = data_train$sample_metadata[positions, ]
    rm(positions)

### Consider only baseline observations for the validation set

    positions = grep("G1", colnames(data_test$feature_table))
    data_test$feature_table = data_test$feature_table[, positions]
    data_test$sample_metadata = data_test$sample_metadata[positions, ]
    rm(positions)

## BayesCoop Implementation

```         
set.seed(1)
result = bayesCoop(data_train, data_test, family = "gaussian", 
                   ss = c(0.05, 1), group = TRUE,
                   alpha_dirich = 1, maxit = 1, 
                   bbiters = 110, bbburn = 10,
                   abd_thresh = 0, prev_thresh = 0.1,
                   Warning = TRUE, verbose = TRUE)

## EM Coordinate Decent Iterations: 1 
## Computational time: 0.022 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.021 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.033 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.039 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.027 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.023 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.015 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.024 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.013 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.023 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.019 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.025 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.019 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.023 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.021 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.026 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.019 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.035 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.021 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.025 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.016 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.028 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.015 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.021 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.019 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.038 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.023 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.029 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.025 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.024 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.018 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.024 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.021 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.017 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.018 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.015 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.031 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.02 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.032 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.016 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.028 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.029 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.026 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.027 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.029 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.026 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.018 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.024 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.021 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.026 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.028 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.024 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.022 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.026 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.022 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.028 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.021 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.026 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.019 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.02 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.017 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.024 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.023 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.024 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.026 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.029 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.017 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.025 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.017 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.022 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.021 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.022 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.017 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.019 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.029 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.021 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.028 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.019 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.022 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.028 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.033 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.017 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.019 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.016 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.023 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.027 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.032 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.018 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.027 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.036 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.02 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.022 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.028 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.017 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.023 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.018 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.019 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.025 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.03 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.018 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.022 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.021 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.021 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.022 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.027 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.024 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.019 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.02 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.02 minutes 
## EM Coordinate Decent Iterations: 1 
## Computational time: 0.027 minutes

#result$beta_postmed

result$rho_postmed

## [1] 0.0001039843

#result$beta_samples

result$y_pred

## [1]  20.352509  -2.883056 -12.267622   5.012925  -4.011487 -13.139735  -6.224429
## [8]  13.160894

result$mspe

## [1] 488.5006

result$time

## [1] 3.172
```

## Visualization of Posterior Samples
Let's check MCMC convergence of the reciprocal Bayesian LASSO estimator through two visualizations: trace plots and histograms.

```         
######################################
# Visualization of Posterior Samples #
######################################

##############
# Trace Plot #
##############

library(coda)
par(mar=c(2,2,2,2))
plot(mcmc(result$beta_samples[,1:9]),density=FALSE,smooth=TRUE)
```
![](https://github.com/himelmallick/BayesCOOP/blob/master/misc/unnamed-chunk-13-1.png)

```         
#############
# Histogram #
#############

library(psych)
multi.hist(result$beta_samples[,1:9],density=TRUE,main="")
```

![](https://github.com/himelmallick/BayesCOOP/blob/master/misc/unnamed-chunk-14-1.png)
