library(BayesCOOP)
set.seed(17)
n <- 36L
x <- list(a = matrix(rnorm(n * 3), n, 3),
          b = matrix(rnorm(n * 2), n, 2))
colnames(x$a) <- paste0("a", 1:3)
colnames(x$b) <- paste0("b", 1:2)
y <- 10 + 2 * x$a[, 1] - x$b[, 2] + rnorm(n, sd = 0.2)
data <- list(
  feature_table = t(do.call(cbind, x)),
  sample_metadata = data.frame(Y = y),
  feature_metadata = data.frame(featureID = unlist(lapply(x, colnames)),
                                featureType = rep(c("a", "b"), c(3, 2))))
foldid <- rep(1:3, times = c(9, 12, 15))
fit <- BayesCOOP(data, bb = FALSE, verbose = FALSE,
                 control = list(alpha = c(0.5, 1), rho = c(0, 0.5),
                                foldid = foldid, s = "lambda.min", keep = TRUE))
cv <- fit$cv_MAP
stopifnot(nrow(cv$alpha.rho.grid) == 4L,
          fit$alpha_MAP %in% c(0.5, 1), fit$rho_MAP %in% c(0, 0.5),
          identical(cv$foldid, as.integer(foldid)),
          length(fit$beta_MAP) == 5L,
          identical(names(fit$beta_MAP), unname(unlist(lapply(x, colnames)))))
# Reproduce held-out scores independently for every alpha/rho combination.
for (result in cv$combo.results) {
  manual <- matrix(NA_real_, length(result$lambda), 3)
  for (f in 1:3) {
    test <- foldid == f
    fold_fit <- multiview::multiview(
      lapply(x, function(z) z[!test, , drop = FALSE]), y[!test],
      family = gaussian(), alpha = result$alpha, rho = result$rho,
      lambda = result$lambda)
    pred <- predict(fold_fit, newx = lapply(x, function(z) z[test, , drop = FALSE]),
                    s = result$lambda, type = "response")
    manual[, f] <- colMeans((y[test] - pred)^2)
  }
  stopifnot(isTRUE(all.equal(manual, result$cvm)),
            result$lambda.min == result$lambda[which.min(apply(manual, 1, weighted.mean, w = tabulate(foldid)))])
}
stopifnot(isTRUE(all.equal(cv$cvm, colMeans((y - cv$fit.preval)^2))))
# Final model must use all observations at the selected settings.
refit <- multiview::multiview(x, y, family = gaussian(), alpha = fit$alpha_MAP,
                             rho = fit$rho_MAP, lambda = fit$lambda_MAP)
stopifnot(isTRUE(all.equal(as.numeric(coef(refit))[-1], unname(fit$beta_MAP))))
# MAP predictions must not depend on held-out outcomes or held-out scaling.
pred <- predict(fit, data, bb = FALSE)$y_pred
changed <- data
changed$sample_metadata$Y <- y + 100
stopifnot(isTRUE(all.equal(pred, predict(fit, changed, bb = FALSE)$y_pred)),
          isTRUE(all.equal(pred, as.numeric(predict(refit, newx = x)))),
          isTRUE(all.equal(pred, as.numeric(fit$intercept_MAP +
                                           do.call(cbind, x) %*% fit$beta_MAP))))
one <- lapply(data, identity)
one$feature_table <- one$feature_table[, 1, drop = FALSE]
one$sample_metadata <- one$sample_metadata[1, , drop = FALSE]
stopifnot(isTRUE(all.equal(predict(fit, one, bb = FALSE)$y_pred, pred[1])))
# Fixed alpha/rho, automatic lambda, 1se selection, and local reproducible seed.
args <- list(alpha = 0.7, rho = 0.25, nfolds = 3, seed = 42)
rng <- .Random.seed
fixed <- BayesCOOP(data, bb = FALSE, verbose = FALSE, control = args)
stopifnot(identical(rng, .Random.seed), fixed$alpha_MAP == 0.7,
          fixed$rho_MAP == 0.25, fixed$lambda_MAP == fixed$cv_MAP$lambda.1se,
          fixed$cv_MAP$lambda.1se >= fixed$cv_MAP$lambda.min)
again <- BayesCOOP(data, bb = FALSE, verbose = FALSE, control = args)
stopifnot(identical(fixed$cv_MAP$foldid, again$cv_MAP$foldid),
          identical(fixed$beta_MAP, again$beta_MAP))
# Explicit lambda grid and invalid control inputs.
explicit <- BayesCOOP(data, bb = FALSE, verbose = FALSE,
                      control = list(rho = 0.5, lambda = c(1, 0.1, 0.01), foldid = foldid))
stopifnot(explicit$lambda_MAP %in% c(1, 0.1, 0.01))
for (bad in list(list(), list(rho = -1), list(rho = 0, alpha = 2),
                 list(rho = 0, foldid = rep(1, n)), list(rho = 0, nfolds = 1))) {
  stopifnot(inherits(try(BayesCOOP(data, bb = FALSE, verbose = FALSE,
                                  control = bad), silent = TRUE), "try-error"))
}
