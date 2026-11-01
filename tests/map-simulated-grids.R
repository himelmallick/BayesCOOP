# Small end-to-end simulations of the spike-and-slab MAP with CV over rho.
library(BayesCOOP)

simulate_data <- function(n, nviews, seed) {
  set.seed(seed)
  x <- lapply(seq_len(nviews), function(v) {
    z <- matrix(rnorm(n * 3), n, 3)
    colnames(z) <- paste0("view", v, "_", seq_len(3))
    z
  })
  y <- 5 + 2 * x[[1]][, 1] - 1.5 * x[[2]][, 2] + rnorm(n, sd = 0.3)
  list(
    feature_table = t(do.call(cbind, x)),
    sample_metadata = data.frame(Y = y),
    feature_metadata = data.frame(
      featureID = unlist(lapply(x, colnames)),
      featureType = rep(paste0("view", seq_len(nviews)), each = 3)))
}

cases <- list(
  list(name = "two views, alpha = 1", nviews = 2, seed = 101,
       control = list(rho = c(0, 0.5, 1), alpha = 1, nfolds = 5, seed = 123)),
  list(name = "three views, alpha = 0.5", nviews = 3, seed = 202,
       control = list(rho = c(0, 0.25, 0.5, 1), alpha = 0.5, nfolds = 5, seed = 123)))

for (case in cases) {
  train <- simulate_data(60, case$nviews, case$seed)
  test <- simulate_data(30, case$nviews, case$seed + 1)
  fit <- BayesCOOP(train, bb = FALSE, verbose = FALSE, warning = FALSE, control = case$control)
  cv <- fit$cv_MAP
  pred <- predict(fit, test, bb = FALSE)
  mse <- mean((pred$y_pred - pred$y_valid)^2)
  baseline <- mean((mean(train$sample_metadata$Y) - pred$y_valid)^2)
  stopifnot(nrow(cv$cv.table) == length(case$control$rho),
            all(is.finite(cv$cv.table$cvm)), all(cv$cv.table$n_failed_folds == 0),
            length(pred$y_pred) == 30L, all(is.finite(pred$y_pred)),
            fit$rho_MAP %in% case$control$rho,
            fit$rho_MAP == cv$cv.table$rho[which.min(cv$cv.table$cvm)],
            fit$alpha_MAP == case$control$alpha,
            is.finite(mse), mse < baseline)
  cat(sprintf("PASS: %s | rho=%g | test MSE=%.4f (mean baseline=%.4f)\n",
              case$name, fit$rho_MAP, mse, baseline))
}
