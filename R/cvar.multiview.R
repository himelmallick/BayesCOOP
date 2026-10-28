# Gaussian adaptation of cvar.multiview for the MAP-only path.
# BayesCOOP currently supports continuous outcomes only.
cvar.multiview <- function(x.list, y, alpha = 1, rho,
                           lambda = NULL, s = "lambda.1se", nfolds = 10,
                           foldid = NULL, keep = FALSE, seed = NULL,
                           verbose = FALSE) {
  if (!is.null(seed)) {
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit(if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }, add = TRUE)
    set.seed(seed)
  }
  s <- match.arg(s, c("lambda.1se", "lambda.min"))
  validate_grid <- function(z, name, lower, upper = Inf) {
    if (!is.numeric(z) || !length(z) || any(!is.finite(z)) ||
        any(z < lower | z > upper)) {
      stop(name, " must contain finite values between ", lower, " and ", upper, ".")
    }
    sort(unique(z))
  }
  alpha <- validate_grid(alpha, "alpha", 0, 1)
  rho <- validate_grid(rho, "rho", 0)
  if (!is.null(lambda)) {
    lambda <- sort(validate_grid(lambda, "lambda", 0), decreasing = TRUE)
  }
  x.list <- lapply(x.list, as.matrix)
  y <- as.numeric(y)
  N <- length(y)
  if (length(x.list) < 2L || N < 3L ||
      any(vapply(x.list, nrow, integer(1)) != N)) {
    stop("MAP cross-validation requires at least two views and three observations with matching rows.")
  }
  if (any(!is.finite(y)) || any(vapply(x.list, function(x) any(!is.finite(x)), logical(1)))) {
    stop("MAP inputs must contain only finite values.")
  }
  if (is.null(foldid)) {
    if (length(nfolds) != 1L || !is.finite(nfolds) || nfolds < 2 || nfolds != as.integer(nfolds)) {
      stop("nfolds must be an integer of at least 2.")
    }
    foldid <- sample(rep(seq_len(min(nfolds, N)), length.out = N))
  } else {
    if (length(foldid) != N || anyNA(foldid)) stop("foldid must have one non-missing value per observation.")
    foldid <- match(foldid, unique(foldid))
  }
  K <- length(unique(foldid))
  if (K < 2L || any(N - tabulate(foldid, K) < 2L)) {
    stop("CV requires at least two folds and two training observations in every fold.")
  }
  fold_weights <- tabulate(foldid, K)
  mv_fit <- function(x, y, a, r, lambda = NULL) {
    args <- list(x_list = x, y = y, family = stats::gaussian(),
                 alpha = a, rho = r, lambda = lambda)
    out <- tryCatch(do.call(multiview::multiview, args), error = identity)
    if (!inherits(out, "error")) return(out)
    if (grepl("cannot correct step size", conditionMessage(out), fixed = TRUE)) {
      retry <- tryCatch(do.call(multiview::multiview,
                                c(args, list(lambda.min.ratio = 0.05,
                                             nlambda = 50, maxit = 1e5))),
                        error = identity)
      if (!inherits(retry, "error")) return(retry)
    }
    stop(conditionMessage(out), call. = FALSE)
  }
  grid <- expand.grid(alpha = alpha, rho = rho)
  results <- vector("list", nrow(grid))
  scores <- rep(NA_real_, nrow(grid))
  failures <- character()
  for (g in seq_len(nrow(grid))) {
    a <- grid$alpha[g]
    r <- grid$rho[g]
    if (verbose) message(sprintf("Tuning alpha = %.3f, rho = %.3f (%d/%d)", a, r, g, nrow(grid)))
    attempt <- tryCatch({
      lambda_g <- if (is.null(lambda)) mv_fit(x.list, y, a, r)$lambda else lambda
      predmat <- matrix(NA_real_, N, length(lambda_g))
      cvm <- matrix(NA_real_, length(lambda_g), K)
      for (f in seq_len(K)) {
        test <- foldid == f
        fit <- mv_fit(lapply(x.list, function(x) x[!test, , drop = FALSE]),
                      y[!test], a, r, lambda_g)
        preds <- as.matrix(stats::predict(fit,
          newx = lapply(x.list, function(x) x[test, , drop = FALSE]),
          s = lambda_g, type = "response"))
        predmat[test, ] <- preds
        cvm[, f] <- colMeans((y[test] - preds)^2)
      }
      cvmean <- apply(cvm, 1, stats::weighted.mean, w = fold_weights)
      cvse <- vapply(seq_along(lambda_g), function(j) {
        sqrt(stats::weighted.mean((cvm[j, ] - cvmean[j])^2, fold_weights) / (K - 1L))
      }, numeric(1))
      idx <- which.min(cvmean)
      lambda_min <- lambda_g[idx]
      lambda_1se <- max(lambda_g[cvmean <= cvmean[idx] + cvse[idx]])
      lambda_choice <- if (s == "lambda.min") lambda_min else lambda_1se
      list(alpha = a, rho = r, lambda = lambda_g, cvm = cvm,
           cvmean = cvmean, cvse = cvse, lambda.min = lambda_min,
           lambda.1se = lambda_1se, lambda.choice = lambda_choice,
           score = cvmean[which.min(abs(lambda_g - lambda_choice))],
           fit.preval = if (keep) predmat else NULL)
    }, error = identity)
    if (inherits(attempt, "error")) {
      failures <- c(failures, conditionMessage(attempt))
      if (verbose) message("Skipping combination: ", conditionMessage(attempt))
      next
    }
    results[[g]] <- attempt
    scores[g] <- attempt$score
  }
  valid <- which(is.finite(scores))
  if (!length(valid)) stop("All alpha/rho combinations failed: ", paste(unique(failures), collapse = "; "))
  best <- valid[which.min(scores[valid])]
  selected <- results[[best]]
  final_fit <- mv_fit(x.list, y, selected$alpha, selected$rho, selected$lambda.choice)
  list(multiview.fit = final_fit, alpha.choice = selected$alpha,
       rho.choice = selected$rho, lambda.choice = selected$lambda.choice,
       lambda.min = selected$lambda.min, lambda.1se = selected$lambda.1se,
       lambda = selected$lambda, cvm = selected$cvmean, cvsd = selected$cvse,
       alpha.rho.grid = grid, combo.results = results, score.choice = scores,
       foldid = foldid, nfolds = K, s = s,
       fit.preval = selected$fit.preval)
}
