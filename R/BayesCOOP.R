#' @title BayesCOOP
#'
#' @description This function implements the BayesCOOP methodology for supervised multimodal
#' integration. It combines jittered group spike-and-slab LASSO regularization with intermediate
#' fusion to enable integrative learning across multiple data modalities. For uncertainty
#' quantification, BayesCOOP applies the Bayesian bootstrap to generate approximate posterior
#' samples by performing maximum a posteriori (MAP) estimation on jittered and resampled
#' multimodal datasets. Currently, only continuous outcomes are supported.
#'
#' @param data_train a list of feature_table, sample_metadata and feature_metadata from training
#' data (unstandardized). See \emph{IntegratedLearner} for more details.
#' @param family currently supports only Gaussian family. Default: "gaussian".
#' @param ss a length-2 numeric vector giving the spike/sLab scales c(s0, s1) with s0 < s1.
#' Default: c(0.05, 1). Used only when \code{bb = TRUE}.
#' @param group logical. If TRUE, predictors are grouped by modality and given a group
#' spike-and-slab prior. If FALSE, no grouping is used. Default: TRUE.
#' Used only when \code{bb = TRUE}.
#' @param bb logical. If TRUE, run full Bayesian bootstrap inference; if FALSE, run MAP
#' estimation using multiview cross-validation with alpha and rho supplied in
#' \code{control}. Default: TRUE.
#' @param alpha_dirich Dirichlet concentration for Bayesian bootstrap weights. Default: 1.
#' @param a_rho,b_rho shape parameters of the Beta(a_rho, b_rho) prior on \eqn{\rho}.
#' @param bbiters number of Bayesian bootstrap iterations. Default: 1100.
#' @param bbburn number of burn-in iterations discarded from the bootstrap. Default: 100.
#' @param maxit maximum EM iterations in the inner optimizer. Default: 100.
#' @param warning logical. If TRUE, emit non-convergence warnings. Default: TRUE.
#' @param verbose logical. If TRUE, print iteration counts and runtime. Default: TRUE.
#' @param control MAP-only cross-validation options. Required element \code{rho}
#' and optional \code{alpha} (default 1) supply fixed values or candidate grids.
#' Lambda is generated automatically and selected by K-fold CV. Optional elements:
#' \code{lambda} (explicit lambda grid), \code{s} (\code{"lambda.1se"}, default,
#' or \code{"lambda.min"}), \code{nfolds} (default 10, capped at sample size),
#' \code{foldid}, \code{seed}, and \code{keep} (retain out-of-fold predictions).
#' The selected model is refitted on all training observations using the
#' multiview elastic-net penalty; \code{ss} and \code{group} do not apply.
#'
#' @return If \code{bb = TRUE}, a list with:
#' \item{beta_samples}{posterior draws of regression coefficients}
#' \item{beta_postmed}{posterior median coefficients}
#' \item{rho_samples}{posterior draws of consensus penalty \eqn{\rho}}
#' \item{rho_postmed}{posterior median of \eqn{\rho}}
#' \item{errVar_samples}{posterior draws of residual variance}
#' \item{errVar_postmed}{posterior median of \eqn{\sigma^2}}
#' \item{time}{runtime in minutes}
#' \item{feature_names}{list of feature names in the training dataset}
#'
#' If \code{bb = FALSE}, a list with:
#' \item{beta_MAP}{MAP estimate of regression coefficients}
#' \item{rho_MAP}{selected agreement parameter from K-fold cross-validation}
#' \item{alpha_MAP}{selected elastic-net mixing parameter}
#' \item{lambda_MAP}{selected regularization parameter}
#' \item{intercept_MAP}{intercept of the full-data multiview refit}
#' \item{multiview_fit}{full-data fitted multiview model, used for MAP predictions}
#' \item{cv_MAP}{cross-validation scores, folds, and tuning results}
#' \item{time}{runtime in minutes}
#' \item{feature_names}{list of feature names in the training dataset}
#'
#' @importFrom glmnet glmnet
#' @importFrom MCMCpack rdirichlet
#' @importFrom rmutil rlaplace
#' @importFrom stats rnorm rgamma quantile median runif rbeta pbeta qbeta
#' @importFrom truncnorm rtruncnorm
#'
#' @export
BayesCOOP <- function(data_train, family = "gaussian",
                      ss = c(0.05, 1), group = TRUE,
                      bb = TRUE, alpha_dirich = 1, a_rho = 1, b_rho = 1,
                      bbiters = 1100, bbburn = 100, maxit = 100,
                      warning = TRUE, verbose = TRUE, control = list()) {
  # Some basic checks
  if(length(data_train) != 3 | any(names(data_train) != c("feature_table", "sample_metadata", "feature_metadata")) == TRUE)
    stop("data_train must be a list of feature_table, sample_metadata and feature_metadata!")
  if(length(ss) != 2 | ss[1] >= ss[2])
    stop("ss must be a two-dimensional vector c(s0, s1) with s0 < s1!")
  if(family != "gaussian")
    stop("Non-gaussian families are currently not supported!")
  if(bbiters <= bbburn)
    stop("(bbiters - bbburn) must be a positive integer!")
  
  output <- list()
  
  #########################################################################
  ######################### data pre-processing ###########################
  #########################################################################
  
  y_train <- gen_datalist(data_train)$y
  xList_train <- gen_datalist(data_train)$xList
  
  feature_names <- lapply(xList_train, names)
  
  ## Scaling the predictors
  xList_train_std <- lapply(xList_train, function(foo) as.data.frame(scale(foo)))
  
  ## centering the response
  y_train <- y_train - mean(y_train)
  
  # Check if group information is provided or not
  if(group == TRUE)
    group <- feature_names else group <- NULL
  
  #########################################################################
  ##################### BayesCOOP Implementation ##########################
  #########################################################################
  
  if(bb) {
    message("Implementing BayesCOOP...")
    
    ## this is the default implementation, i.e, implement jittered bb with sampling rho and error variance
    ## Initialization ##
    rho <- 0.5
    pis <- unlist(lapply(xList_train_std, ncol))
    beta0 <- runif(sum(pis), min = -5, max = 5)
    picums <- cumsum(pis)
    theta <- vector("list", length(xList_train_std))
    theta[[1]] <- beta0[1:picums[1]]
    for(ii in 2:length(xList_train_std)){
      theta[[ii]] <- beta0[(picums[ii-1]+1):picums[ii]]
    }
    errVar <- 1
    
    cts <- 0
    beta_samples <- matrix(NA, bbiters - bbburn, length(beta0))
    rho_samples <- numeric(length = bbiters - bbburn)
    errVar_samples <- numeric(length = bbiters - bbburn)
    
    start.time <- Sys.time()
    for(its in 1:bbiters) {
      
      ###### Update rho ######
      denom <- 0
      for(i in 1:(length(xList_train_std) -1)) {
        for(j in (i + 1):length(xList_train_std)) {
          denom <- denom + sum(drop(as.matrix(xList_train_std[[i]]) %*% theta[[i]] - as.matrix(xList_train_std[[j]]) %*% theta[[j]]) ^ 2) + 1e-5
        }
      }
      # rho <- (truncnorm::rtruncnorm(1, a = 0, b = 1, mean = 0, sd = sqrt(errVar / denom))) ^ 2
      rho <- draw_rho(rho, D = denom, sigma2 = errVar, a = a_rho, b = b_rho)
      
      ## Augmented data for training set using updated rho
      dataAug_train <- data.Augment(y_train, xList_train_std, rho)
      y_aug_train <- dataAug_train$y_aug; x_aug_train <- as.matrix(dataAug_train$x_aug)
      
      
      ###### Update beta ######
      likelihood_weights <- nrow(x_aug_train) * MCMCpack::rdirichlet(1, rep(alpha_dirich, nrow(x_aug_train)))
      
      
      # draw mu_t
      jitter <- rmutil::rlaplace(n = (ncol(x_aug_train)+1), m = 0, s = ss[1]) ## uses the smaller scale, i.e., the scale parameter of spike density
      
      fit_bmlasso_bb <- bmlasso.weighted(x = x_aug_train, y = y_aug_train, family = family,
                                         maxit = maxit, alpha = 1, ss = ss,
                                         lhood_weights = likelihood_weights, jitter = jitter,
                                         group = group, warning = warning, verbose = verbose)
      
      beta_new <- fit_bmlasso_bb$coefficients[-1]
      
      theta[[1]] <- beta_new[1:picums[1]]
      for(kk in 2:length(xList_train_std)){
        theta[[kk]] <- beta_new[(picums[kk-1]+1):picums[kk]]
      }
      
      ###### Update error variance #######
      rss <- sum(drop(y_aug_train - x_aug_train %*% beta_new) ^ 2)
      errVar <- 1 / rgamma(1, shape = (nrow(x_aug_train) + 1) / 2, rate = (rss + 1) / 2)
      
      # store samples
      if(its > bbburn) {
        cts <- cts + 1
        beta_samples[cts, ] <- beta_new
        rho_samples[cts] <- rho
        errVar_samples[cts] <- errVar
      }
      
      message("iter = ", its)
    }
    stop.time <- Sys.time()
    
    # Estimated beta, rho and errVar
    beta_postmed <- apply(beta_samples, 2, function(foo) quantile(foo, 0.5))
    rho_postmed <- median(rho_samples)
    errVar_postmed <- median(errVar_samples)
    time <- as.numeric(round(difftime(stop.time, start.time, units="min"), 3), units = "mins")
    
    
    output$beta_samples <- beta_samples
    colnames(output$beta_samples) <- unlist(feature_names)
    output$beta_postmed <- beta_postmed
    names(output$beta_postmed) <- unlist(feature_names)
    output$rho_samples <- rho_samples
    output$rho_postmed <- rho_postmed
    output$errVar_samples <- errVar_samples
    output$errVar_postmed <- errVar_postmed
    output$time <- time
    output$feature_names <- feature_names
    
  } else {
    
    if(is.null(control$rho)){
      stop("At least one non-null value of rho is needed for MAP estimation!")
    }
    
    message("Implementing MAP estimation...")
    
    start.time <- Sys.time()
    cv_control <- control
    if (is.null(cv_control$alpha)) cv_control$alpha <- 1
    allowed <- c("alpha", "rho", "lambda", "s", "nfolds", "foldid", "keep", "seed")
    unknown <- setdiff(names(cv_control), allowed)
    if (length(unknown)) stop("Unknown MAP control option(s): ", paste(unknown, collapse = ", "))
    fit <- do.call(cvar.multiview,
                   c(list(x.list = xList_train,
                          y = gen_datalist(data_train)$y,
                          verbose = verbose), cv_control))
    stop.time <- Sys.time()
    
    # Estimated beta and rho
    coefficients <- as.numeric(stats::coef(fit$multiview.fit, s = fit$lambda.choice))
    betahat_MAP <- coefficients[-1L]
    rho_MAP <- fit$rho.choice
    time <- as.numeric(round(difftime(stop.time, start.time, units="min"), 3), units = "mins")
    
    output$beta_MAP <- betahat_MAP
    names(output$beta_MAP) <- unlist(feature_names)
    output$rho_MAP <- rho_MAP
    output$alpha_MAP <- fit$alpha.choice
    output$lambda_MAP <- fit$lambda.choice
    output$intercept_MAP <- coefficients[1L]
    output$multiview_fit <- fit$multiview.fit
    output$cv_MAP <- fit
    output$time <- time
    output$feature_names <- feature_names
    
  }
    
  #########################################################################
  ################################ Return #################################
  #########################################################################
  class(output) <- "BayesCOOP"
  output
}
    