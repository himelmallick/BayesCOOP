#' Internal helper functions for BayesCOOP
#'
#' These helpers prepare multimodal data, filter features, construct augmented
#' design matrices for multi-view fusion, and fit weighted spike-and-slab
#' models used internally by \code{bayesCoop()}.
#'
#' Not exported.
#'
#' @keywords internal
#' @name bhglmHelpers
#'
#' @importFrom dplyr filter lag
#' @importFrom stats nobs rnorm rgamma quantile median
#' @importFrom glmnet glmnet
#' @importFrom MCMCpack rdirichlet
#' @importFrom rmutil rlaplace
#' @importFrom truncnorm rtruncnorm
#' @importFrom methods is
#' @importFrom MultiAssayExperiment MultiAssayExperiment experiments intersectColumns
#' @importFrom SummarizedExperiment assay colData
#' @importFrom S4Vectors DataFrame
NULL
###############################################
############## helper functions ###############
###############################################

# Function to reformat MAE data for cooperative learning
gen_datalist <- function(data, y_col = "Y") {
    if (!methods::is(data, "MultiAssayExperiment")) {
        stop("`data` must be a MultiAssayExperiment object.")
    }

    data <- MultiAssayExperiment::intersectColumns(data)
    exp_list <- MultiAssayExperiment::experiments(data)

    if (length(exp_list) < 2) {
        stop("BayesCOOP requires at least two assays in `experiments(data)`.")
    }

    if (is.null(names(exp_list)) || any(names(exp_list) == "")) {
        stop("All assays in `experiments(data)` must be named.")
    }

    sample_metadata <- as.data.frame(SummarizedExperiment::colData(data))
    if (!(y_col %in% colnames(sample_metadata))) {
        stop("Response column `", y_col, "` is missing in colData(data)!")
    }
    if (is.null(rownames(sample_metadata))) {
        stop("colData(data) must have row names corresponding to sample IDs.")
    }

    sample_ids <- rownames(sample_metadata)
    if (length(sample_ids) == 0) {
        stop("No shared samples remain after intersectColumns(data).")
    }

    y <- as.matrix(sample_metadata[[y_col]])
    rownames(y) <- sample_ids

    xList <- vector("list", length(exp_list))
    names(xList) <- names(exp_list)

    for (i in seq_along(exp_list)) {
        assay_i <- exp_list[[i]]
        if (methods::is(assay_i, "SummarizedExperiment")) {
            assay_mat <- as.matrix(SummarizedExperiment::assay(assay_i))
        } else {
            assay_mat <- as.matrix(assay_i)
        }

        if (is.null(colnames(assay_mat))) {
            stop("Each assay in `experiments(data)` must have sample IDs as column names.")
        }
        if (!all(sample_ids %in% colnames(assay_mat))) {
            stop(
                "Assay `", names(exp_list)[i],
                "` is missing at least one sample present in colData(data)."
            )
        }
        assay_mat <- assay_mat[, sample_ids, drop = FALSE]
        xList[[i]] <- as.data.frame(t(assay_mat))
    }

    list(xList = xList, y = y)
}

# Function to reconstruct a MultiAssayExperiment object from (xList, y)
reconstruct_data <- function(xList, y, y_col = "Y") {
  if (length(xList) < 2) {
    stop("xList must contain at least two assays.")
  }

  if (is.null(names(xList)) || any(names(xList) == "")) {
    names(xList) <- paste0("view", seq_along(xList))
  }

  sample_ids <- rownames(xList[[1]])
  if (is.null(sample_ids)) {
    sample_ids <- paste0("sample", seq_len(nrow(xList[[1]])))
  }

  for (i in seq_along(xList)) {
    if (is.null(rownames(xList[[i]]))) {
      rownames(xList[[i]]) <- sample_ids
    }
    if (!identical(rownames(xList[[i]]), sample_ids)) {
      stop("All assays in xList must have the same sample IDs in the same order.")
    }
  }

  assays <- lapply(seq_along(xList), function(i) {
    assay_i <- as.matrix(xList[[i]])
    if (is.null(colnames(assay_i))) {
      colnames(assay_i) <- paste0(names(xList)[i], "_feature", seq_len(ncol(assay_i)))
    }
    if (is.null(rownames(assay_i))) {
      rownames(assay_i) <- sample_ids
    }
    t(assay_i)
  })
  names(assays) <- names(xList)

  if (nrow(y) != length(sample_ids)) {
    stop("`y` must have one value per sample in xList.")
  }

  col_data <- S4Vectors::DataFrame(row.names = sample_ids)
  col_data[[y_col]] <- as.vector(y)

  MultiAssayExperiment::MultiAssayExperiment(
    experiments = assays,
    colData = col_data
  )
}

# Function to have augmented data for cooperative learning application
make_row  <- function(x_list, p_x, pair, rho) {
    sqrt_rho <- sqrt(rho); i <- pair[1]; j <- pair[2];
    result  <- lapply(p_x, function(p) matrix(0, nrow = nrow(x_list[[1L]]), ncol = p))
    result[[i]] <- -sqrt_rho * x_list[[i]]
    result[[j]] <- sqrt_rho * x_list[[j]]
    do.call(cbind, result)
}

data.Augment.2views <- function(y, dataList, rho){
    tdataList <- list()
    tdataList[[1]] <- rbind.data.frame(dataList[[1]], -sqrt(rho) * dataList[[1]])
    tdataList[[2]] <- rbind.data.frame(dataList[[2]], sqrt(rho) * dataList[[2]])
    x_aug <- do.call("cbind", tdataList)
    if(is.null(y)) {
        return(list(x_aug = x_aug))
    } else {
        y_aug <- c(y, rep(0, length(y)))
        return(list(y_aug = y_aug, x_aug = x_aug))
    }
}


data.Augment.mviews <- function(y, dataList, rho){
    p_x <- lapply(dataList, ncol)
    nviews <- length(dataList)
    pairs <- apply(utils::combn(nviews, 2), 2, identity, simplify = FALSE)
    npairs <- length(pairs)
    x <- do.call(cbind, dataList)

    rows <- lapply(pairs, make_row, x_list = dataList, p_x = p_x, rho = rho)
    # matching the column names
    for(i in seq_len(length(rows))){
        colnames(rows[[i]]) <- colnames(x)
    }
    xt <- do.call(rbind, c(list(x), rows))
    
    if(is.null(y)) {
        return(list(x_aug = xt))
    } else {
        yt <- c(y, rep(0, length(pairs) * length(y)))
        return(list(y_aug = yt, x_aug = xt))
    }
    
}

data.Augment <- function(y, dataList, rho){
    if(length(dataList) < 2) {
        stop("BayesCOOP requires at least two modalities.")
    } else if(length(dataList) == 2) {
        dataAug <- data.Augment.2views(y = y, dataList = dataList, rho = rho)
    } else {
        dataAug <- data.Augment.mviews(y = y, dataList = dataList, rho = rho)
    }
    return(dataAug)
}

# Function to filter features according to a fixed abundance and prevalence threshold
filter_features <- function(x, abd_threshold = 0, prev_threshold = 0.1) {
    # x is a data frame
    x <- x[, colMeans(x > abd_threshold) > prev_threshold]
    nzv_x <- caret::nearZeroVar(x, names = TRUE)
    x <- x[, setdiff(names(x), nzv_x)]
    return(x)
}

## bmlasso.fit with likelihood weights for Bayesian Bootstrap approach
# jitter: a location-shift vector to create a random shift in the centering of the prior of beta.
# lwood_weights: randomly samples weights for re-weighting the likelihood function.
# Other arguments are same as bmlasso.fit() function in BhGLM package.
bmlasso.fit.weighted <- function (x, y, family = "gaussian", offset = NULL, epsilon = 1e-04,
                                  maxit = 50, init = rep(0, ncol(x)), alpha = 1,
                                  ss = c(0.05, 1), b = 1, lhood_weights = lhood_weights,
                                  jitter = 0, group = NULL, theta.weights = NULL, inter.hierarchy = NULL,
                                  inter.parents = NULL, warning = FALSE)
{
    ss <- sort(ss)
    ss <- ifelse(ss <= 0, 0.001, ss)
    prior.scale <- ss[length(ss)]
    if (family == "cox") {
        intercept <- FALSE
    } else {intercept <- TRUE}
    x0 <- x
    if (intercept)
        x0 <- cbind(1, x)

    d <- bhglm.prepare(x = x0, intercept = intercept, prior.mean = jitter,
                 prior.sd = 1, prior.scale = prior.scale, prior.df = 1,
                 group = group)
    x <- d$x
    prior.scale <- d$prior.scale

    ###############################################
    prior.mean <- d$prior.mean[-length(d$prior.mean)] # Additional line
    ###############################################

    group <- d$group
    group.vars <- d$group.vars
    ungroup.vars <- d$ungroup.vars
    if (intercept) {
        x <- x[, -1]
        prior.scale <- prior.scale[-1]
    }
    if (length(ss) != 2)
        stop("ss should have two positive values")
    gvars <- unlist(group.vars)
    theta <- p <- rep(0.5, length(gvars))
    names(theta) <- names(p) <- gvars
    if (is.null(theta.weights))
        theta.weights <- rep(1, length(gvars))
    if (length(theta.weights) != length(gvars))
        stop("all grouped variables should have theta.weights")
    if (any(theta.weights > 1 | theta.weights < 0))
        stop("theta.weights should be in [0,1]")
    names(theta.weights) <- gvars
    if (length(b) < length(group.vars))
        b <- c(b, rep(b[length(b)], length(group.vars) - length(b)))
    b <- b[seq_len(length(group.vars))]
    bb <- b
    if (is.null(init)) {
        for (k in seq_len(5)) {
            ps <- ss[1] + (k - 1) * 0.01
            if (family == "cox")
                ps <- min(ss[1] + (k - 1) * 0.01, 0.08)
            alpha0 <- ifelse(alpha == 1, 0.95, 0.05)
            f <- glmnet(x = x, y = y, family = family, offset = offset,
                        alpha = alpha0, lambda = 1/(nrow(x) * ps),
                        weights = lhood_weights, standardize = FALSE)
            b <- as.numeric(f$beta)
            if (any(b != 0))
                break
        }
    } else
        b <- as.numeric(init)
    names(b) <- gvars # names(b) <- colnames(x)
    b <- ifelse(b == 0, 0.001, b)
    init <- b
    devold <- 0
    conv <- FALSE
    for (iter in seq_len(maxit)) {
        if (alpha == 1)
            out <- bhglm.update.scale.p(prior = "mde", b0 = (b[gvars] - prior.mean[-1]),
                                  ss = ss, theta = theta) else
                                      out <- bhglm.update.scale.p(prior = "mt", df = 1e+10,
                                   b0 = b[gvars], ss = ss, theta = theta)
        # prior.scale[gvars] <- out[[1]]
        prior.scale <- out[[1]]
        p <- out[[2]]
        if (!is.matrix(group))
            theta <- bhglm.update.ptheta.group(group.vars = group.vars,
                                         p = p, w = theta.weights, b = bb) else
                                             theta <- bhglm.update.ptheta.network(theta = theta, p = p,
                                            w = group)
        if (!is.null(inter.hierarchy))
            # theta.weights <- bhglm.update.theta.weights(gvars = gvars,
            #                                       theta.weights = theta.weights, inter.hierarchy = inter.hierarchy,
            #                                       inter.parents = inter.parents, p = p)
            theta.weights <- bhglm.update.theta.weights(eff.hierarchy = inter.hierarchy,
                eff.parents = inter.parents, p = p)
        Pf <- 1/(prior.scale + 1e-10)
        f <- glmnet(x = x, y = y, family = family, offset = offset,
                    alpha = alpha, penalty.factor = Pf, weights = lhood_weights,
                    lambda = sum(Pf)/(nrow(x) * ncol(x)), standardize = FALSE)
        b <- as.numeric(f$beta)
        names(b) <- gvars                   # names(b) <- colnames(x)
        dev <- deviance(f)
        if (abs(dev - devold)/(0.1 + abs(dev)) < epsilon & iter >
            5) {
            conv <- TRUE
            break
        } else devold <- dev
    }
    if (warning & !conv)
        warning("algorithm did not converge", call. = FALSE)
    f$x <- x
    f$y <- y
    f$family <- family
    f$ss <- ss
    f$coefficients <- as.numeric(coef(f))
    names(f$coefficients) <- rownames(coef(f))
    f$linear.predictors <- predict(f, newx = x, type = "link",
                                   offset = offset)
    # if (family == "gaussian")
    #     f$dispersion <- bhglm.bglm(y ~ f$linear.predictors - 1, start = 1,
    #                          prior = De(1, 0), verbose = FALSE)$dispersion
    f$iter <- iter
    f$init <- init
    f$aic <- deviance(f) + 2 * f$df
    f$offset <- offset
    f$prior.scale <- prior.scale
    f$penalty.factor <- Pf
    f$group <- group
    f$group.vars <- group.vars
    f$ungroup.vars <- ungroup.vars
    f$p <- p
    f$ptheta <- theta
    f$b <- bb
    f$theta.weights <- theta.weights
    return(f)
}

## bmlasso with likelihood weights for Bayesian Bootstrap approach
# jitter: a location-shift vector to create a random shift in the centering of the prior of beta.
# lwood_weights: randomly samples weights for re-weighting the likelihood function.
# Other arguments are same as bmlasso() function in BhGLM package.
bmlasso.weighted <- function (x, y, family = c("gaussian", "binomial", "poisson",
                           "cox"), offset = NULL, epsilon = 1e-04, maxit = 50, init = NULL,
          alpha = c(1, 0), ss = c(0.04, 0.5), b = 1, group = NULL, lhood_weights = lhood_weights, jitter = NULL,
          theta.weights = NULL, inter.hierarchy = NULL, inter.parents = NULL,
          warning = FALSE, verbose = FALSE)
{
    start.time <- Sys.time()
    call <- match.call()
    x <- as.matrix(x)
    if (is.null(colnames(x)))
        colnames(x) <- paste("x", seq_len(ncol(x)), sep = "")
    nobs <- nrow(x)
    if (NROW(y) != nobs)
        stop("nobs of 'x' and 'y' are different")
    inc <- apply(cbind(y, x), 1, function(z) !any(is.na(z)))
    if (!is.null(offset)) {
        if (length(offset) != nobs)
            stop("nobs of 'x' and 'offset' are different")
        inc <- apply(cbind(y, x, offset), 1, function(z) !any(is.na(z)))
    }
    y <- y[inc]
    x <- x[inc, ]
    offset <- offset[inc]
    family <- family[1]
    if (family == "cox")
        if (!is.Surv(y))
            stop("'y' should be a 'Surv' object")
    if (!is.null(init) & length(init) != ncol(x))
        stop("give an initial value to each coefficient (not intercept)")
    alpha <- alpha[1]
    f <- bmlasso.fit.weighted(x = x, y = y, family = family, offset = offset,
                     epsilon = epsilon, maxit = maxit, init = init, group = group,
                     alpha = alpha, ss = ss, b = b, theta.weights = theta.weights,
                     lhood_weights = lhood_weights, jitter = jitter,
                     inter.hierarchy = inter.hierarchy, inter.parents = inter.parents,
                     warning = warning)
    f$call <- call
    if (family == "cox")
        class(f) <- c(class(f), "bmlasso", "COXPH")
    else class(f) <- c(class(f), "bmlasso", "GLM")
    stop.time <- Sys.time()
    minutes <- round(difftime(stop.time, start.time, units = "min"),
                     3)
    if (verbose) {
        message("EM Coordinate Decent Iterations: ", f$iter)
        message("Computational time: ", minutes, " minutes")
    }
    return(f)
}

# Function to choose optimal rho for bmlasso
bmlasso_cv <- function(y, xList,
                      rho_grid, group, ss, family,
                      maxit, warning, verbose){

    if (length(xList) < 2) {
        stop("BayesCOOP requires at least two modalities.")
    }

    train_indx <- sample(seq_len(length(y)), round(length(y) * 0.7))
    test_indx <- setdiff(seq_len(length(y)), train_indx)
    y_train <- y[train_indx]; y_train <- y_train - mean(y_train)
    y_test <- y[test_indx]; y_test <- y_test - mean(y_test)
    xList_train <- lapply(xList, function(foo) foo[train_indx, ])
    xList_test <- lapply(xList, function(foo) foo[test_indx, ])

    ## Scaling the predictors
    xList_train_std <- lapply(xList_train, function(foo) as.data.frame(scale(foo)))
    xList_test_std <- lapply(xList_test, function(foo) as.data.frame(scale(foo)))

    # initialize variables to store output
    mspe_bmlasso_group <- rep(0, length(rho_grid))
    beta_hat_group <- vector("list", length(rho_grid))

    # grid search for optimal rho
    for(i in seq_len(length(rho_grid))){
        ## Augmented data
        if(length(xList_train) == 2) {
            dataAug_train <- data.Augment.2views(y_train, xList_train_std, rho_grid[i])
            dataAug_test <- data.Augment.2views(y_test, xList_test_std, rho_grid[i])
        } else {
            dataAug_train <- data.Augment.mviews(y_train, xList_train_std, rho_grid[i])
            dataAug_test <- data.Augment.mviews(y_test, xList_test_std, rho_grid[i])
        }

        ## bmlasso with grouped predictors for fixed rho
        y_aug_train <- dataAug_train$y_aug; x_aug_train <- dataAug_train$x_aug
        fit_bmlasso_group <- bhglm.bmlasso(x = as.matrix(x_aug_train), y = y_aug_train, family = family,
                                    maxit = maxit, alpha = 1, ss = ss,
                                    group = group, warning = TRUE, verbose = TRUE)
        beta_hat_group[[i]] <- fit_bmlasso_group$coefficients[-1]

        ## Prediction
        y_aug_test <- dataAug_test$y_aug; x_aug_test <- dataAug_test$x_aug
        y_aug_hat_group <- as.matrix(x_aug_test) %*% beta_hat_group[[i]]

        ## Consider the first length(test.y) values from yhat[[k]] for predicted values of y
        y_pred_group <- y_aug_hat_group[seq_len(length(y_test))]

        ## mean squared prediction error calculation
        mspe_bmlasso_group[i] <- mean((y_pred_group - y_test) ^ 2)
    }
    # find index corresponding to lowest mse
    opt_index <- which.min(mspe_bmlasso_group)
    rho_bmlasso <- rho_grid[opt_index]

    ### fitting optimum model to obtain beta_hat
    y <- y - mean(y)
    xList_std <- lapply(xList, function(foo) as.data.frame(scale(foo)))

    if(length(xList_std) == 2) {
        dataAug_train <- data.Augment.2views(y, xList_std, rho_bmlasso)
    } else {
        dataAug_train <- data.Augment.mviews(y, xList_std, rho_bmlasso)
    }

    yAug <- dataAug_train$y_aug; xAug <- dataAug_train$x_aug
    model_bmlasso <- bhglm.bmlasso(x = as.matrix(xAug), y = yAug, family = family,
                            maxit = maxit, alpha = 1, ss = ss,
                            group = group, warning = TRUE, verbose = TRUE)
    beta_hat_bmlasso <- model_bmlasso$coefficients[-1]

    return(list(beta_hat_bmlasso = beta_hat_bmlasso, rho_bmlasso = rho_bmlasso))
}
