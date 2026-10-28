# Small end-to-end simulations for explicit and automatic lambda grids.
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
  list(name = "two views, explicit lambda grid", nviews = 2, seed = 101,
       control = list(rho = c(0, 0.5, 1), alpha = c(0.5, 0.75, 1),
                      lambda = c(1, 0.1, 0.01), nfolds = 5,
                      s = "lambda.min", seed = 123)),
  list(name = "three views, automatic lambda grid", nviews = 3, seed = 202,
       control = list(rho = c(0, 0.5, 1), alpha = c(0.5, 0.75, 1),
                      nfolds = 5, s = "lambda.1se", seed = 123)))

for (case in cases) {
  train <- simulate_data(60, case$nviews, case$seed)
  test <- simulate_data(30, case$nviews, case$seed + 1)
  fit <- BayesCOOP(train, bb = FALSE, verbose = FALSE, control = case$control)
  cv <- fit$cv_MAP
  pred <- predict(fit, test, bb = FALSE)
  mse <- mean((pred$y_pred - pred$y_valid)^2)
  baseline <- mean((mean(train$sample_metadata$Y) - pred$y_valid)^2)
  stopifnot(length(cv$combo.results) == 9L,
            all(vapply(cv$combo.results, function(z) !is.null(z), logical(1))),
            length(pred$y_pred) == 30L, all(is.finite(pred$y_pred)),
            fit$rho_MAP %in% case$control$rho,
            fit$alpha_MAP %in% case$control$alpha,
            is.finite(mse), mse < baseline)
  if (!is.null(case$control$lambda)) {
    stopifnot(fit$lambda_MAP %in% case$control$lambda)
  } else {
    stopifnot(all(vapply(cv$combo.results, function(z) length(z$lambda) > 1L, logical(1))))
  }
  best <- which.min(cv$score.choice)
  chosen <- cv$combo.results[[best]]
  stopifnot(fit$rho_MAP == chosen$rho, fit$alpha_MAP == chosen$alpha,
            fit$lambda_MAP == chosen$lambda.choice)
  cat(sprintf("PASS: %s | alpha=%g rho=%g lambda=%.6g | test MSE=%.4f (mean baseline=%.4f)\n",
              case$name, fit$alpha_MAP, fit$rho_MAP, fit$lambda_MAP, mse, baseline))
}
