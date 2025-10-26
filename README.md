# BayesCOOP: Bayesian Cooperative Learning for Multimodal Integration

---

**BayesCOOP** is an R package for supervised multimodal integration with **continuous outcomes**.  

It provides a scalable Bayesian alternative to traditional multiview regression by:  
- enforcing group-wise sparsity via a spike-and-slab Laplace prior to identify informative features and modalities, and  
- quantifying predictive uncertainty through the Bayesian bootstrap.  

At its core, BayesCOOP combines fast maximum a posteriori (MAP) estimation with Bayesian bootstrap–based uncertainty quantification under a *jittered* group spike-and-slab double exponential prior, enabling simultaneous prediction and feature selection with principled uncertainty estimates using pseudo-posteriors. 

The package supports two modes:  
1. **Full Bayesian mode (`bb = TRUE`)** — generates posterior draws, credible intervals, and an agreement parameter reflecting cross-view cooperation.  
2. **Fast MAP mode (`bb = FALSE`)** — provides cross-validated MAP estimates without full posterior sampling for improved speed.  

---

## ⚙️ Installation

You can install the development version directly from GitHub:

```r
install.packages("devtools")
devtools::install_github("himelmallick/BayesCOOP")
library(BayesCOOP)
```

---

## 🚀 Quick Example

Here is a minimal call to `BayesCOOP()`:

```r
fit <- BayesCOOP::BayesCOOP(
  data_train,           # list with feature_table, sample_metadata and feature_metadata (train)
  data_test,            # same structure for held-out / validation data
  family       = "gaussian",   # currently only continuous outcomes are supported
  ss           = c(0.05, 1),   # spike-and-slab scales (s0, s1); controls sparsity
  group        = TRUE,         # whether to apply group-wise shrinkage by modality
  bb           = TRUE,         # TRUE = Bayesian bootstrap mode (default); FALSE = MAP only
  alpha_dirich = 1,            # Dirichlet weight parameter for the Bayesian bootstrap
  bbiters      = 1100,         # total bootstrap iterations
  bbburn       = 100,          # burn-in iterations discarded from the bootstrap
  maxit        = 100,          # max EM iterations for the inner optimizer
  filter       = TRUE,         # filter low-abundance / low-prevalence features first
  abd_thresh   = 0,            # abundance threshold for filtering
  prev_thresh  = 0.1,          # prevalence threshold for filtering
  Warning      = TRUE,         # print convergence warnings
  verbose      = TRUE,         # print iteration progress and timing
  control      = list()        # (used if bb = FALSE; see below)
)
```

### What do I get back?

When `bb = TRUE` (the default Bayesian mode), the returned object includes:

- `mspe`: mean squared prediction error on the test set  
- `y_pred`: predicted response values  
- `y_samples`: posterior predictive draws  
- `beta_postmed`: posterior median of the coefficients  
- `rho_postmed`: posterior median of the agreement parameter across views  
- `time`: runtime in minutes

When `bb = FALSE`, BayesCOOP runs fast MAP estimation instead of the Bayesian bootstrap. In that case you must pass `control = list(rho = c(...))`, e.g.:

```r
fit_map <- BayesCOOP::BayesCOOP(
  data_train,
  data_test,
  bb      = FALSE,
  control = list(rho = c(0, 0.5, 1))  # candidate cross-view agreement values
)

fit_map$mspe     # MSPE for the best fusion level
fit_map$rho_MAP  # selected agreement parameter
fit_map$beta_MAP # MAP coefficients
```

In both modes, `mspe` is the main out-of-sample accuracy metric you can report.

---

## 📘 Full Tutorial

For an in-depth workflow — including real data preprocessing, baseline comparisons, and performance benchmarking — please see the full tutorial:

📄 [View the BayesCOOP Tutorial](https://raw.githack.com/himelmallick/BayesCOOP/master/vignettes/BayesCOOP.html)

---

## 📚 Citation

If you use **BayesCOOP** in your work, please cite:

> Roy, S., Sarkar, S., Paul, E., Basak, P., Yi, N., & Mallick, H. (2025).  
> **Bayesian Cooperative Learning for Multimodal Integration.** *bioRxiv.*  
> [https://doi.org/10.1101/2025.10.23.684056](https://doi.org/10.1101/2025.10.23.684056)

You can also use the standard R citation command:

```r
citation("BayesCOOP")
```

---

## 🐞 Issues

We are happy to troubleshoot any issues with the package. Please contact the authors via email or open an issue in the [GitHub repository](https://github.com/himelmallick/BayesCOOP/issues).
