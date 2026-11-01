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
rho_grid <- c(0, 0.5)
fit <- BayesCOOP(data, bb = FALSE, verbose = FALSE, warning = FALSE,
                 control = list(alpha = 0.5, rho = rho_grid, foldid = foldid))
cv <- fit$cv_MAP
stopifnot(nrow(cv$cv.table) == 2L, fit$rho_MAP %in% rho_grid, fit$alpha_MAP == 0.5,
          identical(cv$foldid, as.integer(foldid)), cv$nfolds == 3L,
          length(fit$beta_MAP) == 5L,
          identical(names(fit$beta_MAP), unname(unlist(lapply(x, colnames)))))

# Reproduce every held-out score with the package's own building blocks:
# standardize on the fold's training rows, augment, fit bmlasso, predict.
bmlasso <- getFromNamespace("bhglm.bmlasso", "BayesCOOP")
augment <- getFromNamespace("data.Augment", "BayesCOOP")
fit_manual <- function(xs, yy, rho) {
  mu <- lapply(xs, colMeans)
  sdv <- lapply(xs, function(z) apply(z, 2, sd))
  xstd <- Map(function(z, m, s) as.data.frame(sweep(sweep(z, 2, m), 2, s, "/")), xs, mu, sdv)
  aug <- augment(yy - mean(yy), xstd, rho)
  xa <- data.matrix(aug$x_aug)
  colnames(xa) <- unlist(lapply(xs, colnames))
  f <- bmlasso(x = xa, y = as.numeric(aug$y_aug), family = "gaussian", maxit = 100,
               alpha = 0.5, ss = c(0.05, 1), group = unname(lapply(xs, colnames)),
               warning = FALSE, verbose = FALSE)
  list(beta = as.numeric(f$coefficients[-1]), mu = mu, sdv = sdv, ym = mean(yy))
}
pred_manual <- function(f, xs) {
  as.numeric(do.call(cbind, Map(function(z, m, s) sweep(sweep(z, 2, m), 2, s, "/"),
                                xs, f$mu, f$sdv)) %*% f$beta) + f$ym
}
for (i in seq_along(rho_grid)) {
  manual <- vapply(1:3, function(k) {
    test <- foldid == k
    f <- fit_manual(lapply(x, function(z) z[!test, , drop = FALSE]), y[!test], rho_grid[i])
    mean((y[test] - pred_manual(f, lapply(x, function(z) z[test, , drop = FALSE])))^2)
  }, numeric(1))
  stopifnot(isTRUE(all.equal(manual, unname(cv$cvm.folds[i, ]))),
            isTRUE(all.equal(weighted.mean(manual, tabulate(foldid)), cv$cv.table$cvm[i])))
}
stopifnot(fit$rho_MAP == cv$cv.table$rho[which.min(cv$cv.table$cvm)])

# Final model must use all observations at the selected rho.
refit <- fit_manual(x, y, fit$rho_MAP)
stopifnot(isTRUE(all.equal(refit$beta, unname(fit$beta_MAP))))

# MAP predictions must not depend on held-out outcomes or held-out scaling.
pred <- predict(fit, data, bb = FALSE)$y_pred
changed <- data
changed$sample_metadata$Y <- y + 100
stopifnot(isTRUE(all.equal(pred, predict(fit, changed, bb = FALSE)$y_pred)),
          isTRUE(all.equal(pred, pred_manual(refit, x))))
one <- lapply(data, identity)
one$feature_table <- one$feature_table[, 1, drop = FALSE]
one$sample_metadata <- one$sample_metadata[1, , drop = FALSE]
stopifnot(isTRUE(all.equal(predict(fit, one, bb = FALSE)$y_pred, pred[1])))

# Local reproducible seed: same folds and fit, global RNG left untouched.
args <- list(alpha = 1, rho = 0.25, nfolds = 3, seed = 42)
rng <- .Random.seed
fixed <- BayesCOOP(data, bb = FALSE, verbose = FALSE, warning = FALSE, control = args)
stopifnot(identical(rng, .Random.seed), fixed$alpha_MAP == 1, fixed$rho_MAP == 0.25)
again <- BayesCOOP(data, bb = FALSE, verbose = FALSE, warning = FALSE, control = args)
stopifnot(identical(fixed$cv_MAP$foldid, again$cv_MAP$foldid),
          identical(fixed$beta_MAP, again$beta_MAP))

# Invalid control inputs, including the removed lambda/s options.
for (bad in list(list(), list(rho = -1), list(rho = 0, alpha = 2),
                 list(rho = 0, alpha = c(0.5, 1)),
                 list(rho = 0, foldid = rep(1, n)), list(rho = 0, nfolds = 1),
                 list(rho = 0, lambda = 0.1), list(rho = 0, s = "lambda.min"))) {
  stopifnot(inherits(try(BayesCOOP(data, bb = FALSE, verbose = FALSE,
                                  control = bad), silent = TRUE), "try-error"))
}
cat("map-cv: all checks passed\n")
