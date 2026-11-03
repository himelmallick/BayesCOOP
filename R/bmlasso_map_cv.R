# MAP estimation for BayesCOOP (bb = FALSE).
#
# The model is BayesCOOP's own: a (group) spike-and-slab prior fitted by
# bhglm.bmlasso() on the rho-augmented data. The agreement parameter rho is
# chosen by K-fold cross-validation, organised like cvar.multiview(): every
# candidate rho is fitted on K - 1 folds and scored on the held-out fold, the
# fold errors are averaged (weighted by fold size), the rho with the lowest
# mean CV error is selected, and the model is refitted on all observations at
# that rho. There is no lambda to tune: the penalty inside bmlasso comes from
# the spike/slab scales ss.
#
# Predictors are centred and scaled, and the response centred, using only the
# observations being fitted; held-out folds are transformed with those values.
#
# method = "holdout" instead reproduces the original BayesCOOP MAP tuning
# (the former bmlasso_cv()): one random split of the training data into
# train_frac for fitting and the rest for validation, rho chosen by the
# validation MSPE, then a refit on all observations. As in the original, the
# validation predictors are standardized and the validation response centred
# with their own means/SDs.

bmlasso_map_cv <- function(x.list, y, rho, alpha = 1, nfolds = 10,
                           foldid = NULL, seed = NULL,
                           method = c("kfold", "holdout"), train_frac = 0.7,
                           ss = c(0.05, 1),
                           group = TRUE, family = "gaussian", maxit = 100,
                           warning = FALSE, verbose = FALSE) {
  method <- match.arg(method)
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
  if (!is.numeric(rho) || !length(rho) || any(!is.finite(rho)) || any(rho < 0)) {
    stop("rho must contain finite, non-negative values.")
  }
  rho <- sort(unique(rho))
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) ||
      alpha < 0 || alpha > 1) {
    stop("alpha must be a single value between 0 and 1.")
  }

  x.list <- lapply(x.list, function(x) {
    x <- as.matrix(x)
    storage.mode(x) <- "double"
    x
  })
  y <- as.numeric(y)
  N <- length(y)
  if (length(x.list) < 2L || N < 3L ||
      any(vapply(x.list, nrow, integer(1)) != N)) {
    stop("MAP cross-validation requires at least two views and three observations with matching rows.")
  }
  if (any(!is.finite(y)) || any(vapply(x.list, function(x) any(!is.finite(x)), logical(1)))) {
    stop("MAP inputs must contain only finite values.")
  }
  view_cols <- lapply(seq_along(x.list), function(v) {
    cn <- colnames(x.list[[v]])
    if (is.null(cn)) cn <- paste0("view", v, "_x", seq_len(ncol(x.list[[v]])))
    cn
  })
  names(view_cols) <- names(x.list)
  for (v in seq_along(x.list)) colnames(x.list[[v]]) <- view_cols[[v]]
  groups <- if (isTRUE(group)) unname(view_cols) else NULL

  fit_scaler <- function(x) {
    sig <- apply(x, 2L, stats::sd)
    sig[!is.finite(sig) | sig <= 0] <- 1
    list(center = colMeans(x), scale = sig)
  }
  apply_scaler <- function(x, sc) {
    sweep(sweep(x, 2L, sc$center, "-"), 2L, sc$scale, "/")
  }

  # Fit the augmented spike-and-slab model at one rho.
  fit_one <- function(xs, yy, r) {
    scalers <- lapply(xs, fit_scaler)
    x_std <- Map(function(x, sc) as.data.frame(apply_scaler(x, sc), check.names = FALSE),
                 xs, scalers)
    y_mean <- mean(yy)
    aug <- data.Augment(yy - y_mean, x_std, r)
    x_aug <- data.matrix(aug$x_aug)
    colnames(x_aug) <- unlist(view_cols, use.names = FALSE)
    fit <- bhglm.bmlasso(x = x_aug, y = as.numeric(aug$y_aug), family = family,
                         maxit = maxit, alpha = alpha, ss = ss, group = groups,
                         warning = warning, verbose = FALSE)
    beta <- as.numeric(fit$coefficients[-1L])
    names(beta) <- unlist(view_cols, use.names = FALSE)
    list(beta = beta, scalers = scalers, y_mean = y_mean)
  }
  predict_one <- function(f, xs) {
    x_new <- do.call(cbind, Map(apply_scaler, xs, f$scalers))
    as.numeric(x_new %*% f$beta) + f$y_mean
  }
  final_output <- function(best, cv_table, extra) {
    final <- fit_one(x.list, y, rho[best])
    c(list(beta = final$beta, rho.choice = rho[best], alpha = alpha, method = method,
           x.center = lapply(final$scalers, `[[`, "center"),
           x.scale = lapply(final$scalers, `[[`, "scale"),
           y.mean = final$y_mean, cv.table = cv_table), extra)
  }

  if (method == "holdout") {
    if (length(train_frac) != 1L || !is.finite(train_frac) ||
        train_frac <= 0 || train_frac >= 1) {
      stop("train_frac must be a single value strictly between 0 and 1.")
    }
    train_index <- sample(seq_len(N), round(N * train_frac))
    valid <- setdiff(seq_len(N), train_index)
    if (length(train_index) < 2L || length(valid) < 2L) {
      stop("The holdout split needs at least two observations on each side.")
    }
    x_valid <- lapply(x.list, function(x) x[valid, , drop = FALSE])
    x_valid_std <- do.call(cbind, lapply(x_valid, function(x) apply_scaler(x, fit_scaler(x))))
    y_valid <- y[valid] - mean(y[valid])
    mspe <- rep(NA_real_, length(rho))
    failures <- character()
    for (i in seq_along(rho)) {
      if (verbose) message(sprintf("Tuning rho = %.3f (%d/%d)", rho[i], i, length(rho)))
      mspe[i] <- tryCatch({
        f <- fit_one(lapply(x.list, function(x) x[train_index, , drop = FALSE]),
                     y[train_index], rho[i])
        mean((as.numeric(x_valid_std %*% f$beta) - y_valid)^2)
      }, error = function(e) {
        failures <<- c(failures, conditionMessage(e))
        NA_real_
      })
    }
    if (!any(is.finite(mspe))) {
      stop("All rho values failed on the holdout split: ", paste(unique(failures), collapse = "; "))
    }
    cv_table <- data.frame(rho = rho, cvm = mspe, cvsd = NA_real_,
                           n_failed_folds = as.integer(!is.finite(mspe)))
    return(final_output(which.min(mspe), cv_table,
                        list(train_index = train_index, train_frac = train_frac)))
  }

  # Folds drawn as in cvar.multiview().
  if (is.null(foldid)) {
    if (length(nfolds) != 1L || !is.finite(nfolds) || nfolds < 2 ||
        nfolds != as.integer(nfolds)) {
      stop("nfolds must be an integer of at least 2.")
    }
    foldid <- sample(rep(seq_len(min(nfolds, N)), length.out = N))
  } else {
    if (length(foldid) != N || anyNA(foldid)) {
      stop("foldid must have one non-missing value per observation.")
    }
    foldid <- match(foldid, unique(foldid))
  }
  K <- length(unique(foldid))
  if (K < 2L || any(N - tabulate(foldid, K) < 2L)) {
    stop("CV requires at least two folds and two training observations in every fold.")
  }
  fold_weights <- tabulate(foldid, K)

  cv <- matrix(NA_real_, length(rho), K)
  failures <- character()
  for (i in seq_along(rho)) {
    if (verbose) message(sprintf("Tuning rho = %.3f (%d/%d)", rho[i], i, length(rho)))
    for (k in seq_len(K)) {
      test <- foldid == k
      cv[i, k] <- tryCatch({
        f <- fit_one(lapply(x.list, function(x) x[!test, , drop = FALSE]), y[!test], rho[i])
        pred <- predict_one(f, lapply(x.list, function(x) x[test, , drop = FALSE]))
        mean((y[test] - pred)^2)
      }, error = function(e) {
        failures <<- c(failures, conditionMessage(e))
        NA_real_
      })
    }
  }

  cv_table <- do.call(rbind, lapply(seq_along(rho), function(i) {
    ok <- is.finite(cv[i, ])
    m <- if (any(ok)) stats::weighted.mean(cv[i, ok], fold_weights[ok]) else NA_real_
    se <- if (sum(ok) > 1L) {
      sqrt(stats::weighted.mean((cv[i, ok] - m)^2, fold_weights[ok]) / (sum(ok) - 1L))
    } else NA_real_
    data.frame(rho = rho[i], cvm = m, cvsd = se, n_failed_folds = sum(!ok))
  }))
  if (!any(is.finite(cv_table$cvm))) {
    stop("All rho values failed in cross-validation: ", paste(unique(failures), collapse = "; "))
  }
  final_output(which.min(cv_table$cvm), cv_table,
               list(cvm.folds = cv, foldid = foldid, nfolds = K))
}
