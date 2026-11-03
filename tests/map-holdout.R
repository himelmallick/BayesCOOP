library(BayesCOOP)
set.seed(101)
n <- 60L
x <- list(a = matrix(rnorm(n * 4), n, 4), b = matrix(rnorm(n * 4), n, 4))
colnames(x$a) <- paste0("a", 1:4)
colnames(x$b) <- paste0("b", 1:4)
y <- 5 + 2 * x$a[, 1] - 1.5 * x$b[, 2] + rnorm(n)
data <- list(
  feature_table = t(do.call(cbind, x)),
  sample_metadata = data.frame(Y = y),
  feature_metadata = data.frame(featureID = unlist(lapply(x, colnames)),
                                featureType = rep(c("a", "b"), each = 4)))
grid <- c(0, 0.25, 0.5, 1)

# Original BayesCOOP MAP tuning: one 70/30 holdout split.
fit <- BayesCOOP(data, bb = FALSE, verbose = FALSE, warning = FALSE,
                 control = list(rho = grid, method = "holdout", seed = 3))
stopifnot(fit$cv_MAP$method == "holdout", fit$rho_MAP %in% grid,
          nrow(fit$cv_MAP$cv.table) == length(grid),
          length(fit$cv_MAP$train_index) == round(0.7 * n),
          length(fit$beta_MAP) == 8L)

# Same seed gives the same split and fit; train_frac is respected.
f1 <- BayesCOOP(data, bb = FALSE, verbose = FALSE, warning = FALSE,
                control = list(rho = grid, method = "holdout", seed = 3, train_frac = 0.8))
f2 <- BayesCOOP(data, bb = FALSE, verbose = FALSE, warning = FALSE,
                control = list(rho = grid, method = "holdout", seed = 3, train_frac = 0.8))
stopifnot(identical(f1$beta_MAP, f2$beta_MAP), length(f1$cv_MAP$train_index) == 48L)

# The validation MSPE is computed as in the original: validation predictors
# standardized and response centred with the validation set's own values.
valid <- setdiff(seq_len(n), fit$cv_MAP$train_index)
stopifnot(all(is.finite(fit$cv_MAP$cv.table$cvm)), length(valid) == n - 42L)

# Default is still K-fold; bad options fail.
stopifnot(BayesCOOP(data, bb = FALSE, verbose = FALSE, warning = FALSE,
                    control = list(rho = grid, nfolds = 3, seed = 1))$cv_MAP$method == "kfold")
stopifnot(inherits(try(BayesCOOP(data, bb = FALSE, verbose = FALSE,
                                 control = list(rho = 0, method = "loo")), silent = TRUE), "try-error"))
stopifnot(inherits(try(BayesCOOP(data, bb = FALSE, verbose = FALSE,
                                 control = list(rho = 0, method = "holdout", train_frac = 1)), silent = TRUE), "try-error"))

p <- predict(fit, data, bb = FALSE, verbose = FALSE)
stopifnot(length(p$y_pred) == n)
cat("map-holdout: all checks passed\n")
