#' @title predict.BayesCOOP
#'
#' @description This function predicts the response from a fitted object using BayesCOOP methodology. Currently, only continuous outcomes are supported.
#'
#' @param object fitted object returned by the function \code{\link{BayesCOOP}} implemented on the training dataset.
#' @param newdata a list of feature_table, sample_metadata and feature_metadata from validation dataset (unstandardized). See \emph{IntegratedLearner} for more details.
#' @param family currently supports only Gaussian family. Default: "gaussian".
#' @param bb logical. If TRUE, posterior predictive sampling is performed using the posterior median of \eqn{rho} stored in the fitted object; if FALSE, point prediction is provided for MAP estimate of \eqn{rho} stored in the fitted object. Default: TRUE.
#' @param warning logical. If TRUE, emit non-convergence warnings. Default: TRUE.
#' @param verbose logical. If TRUE, print iteration counts and runtime. Default: TRUE.
#' @param ... further arguments passed from \code{predict()}.
#'
#' @return If \code{bb = TRUE}, a list with:
#' \item{y_samples}{posterior predictive draws (approximate)}
#' \item{y_pred}{posterior predictive median for each held-out sample}
#' \item{y_valid}{True values of response in the validation dataset}
#'
#' If \code{bb = FALSE}, a list with:
#' \item{y_pred}{predicted response values for data}
#' \item{y_valid}{True values of response in the validation dataset}
#'
#' @importFrom stats rnorm quantile median
#' @method predict BayesCOOP
#' @export
predict.BayesCOOP <- function(object, newdata,
                              family = "gaussian",
                              bb = TRUE,
                              warning = TRUE, verbose = TRUE, ...) {
  
  # Some basic checks
  if(length(newdata) != 3 | any(names(newdata) != c("feature_table", "sample_metadata", "feature_metadata")) == TRUE)
    stop("newdata must be a list of feature_table, sample_metadata and feature_metadata!")
  if(family != "gaussian")
    stop("Non-gaussian families are currently not supported!")
  
  output <- list()
  
  #########################################################################
  ######################### data pre-processing ###########################
  #########################################################################
  
  y_test <- gen_datalist(newdata)$y
  xList_test <- gen_datalist(newdata)$xList
  
  ## Scaling the predictors
  xList_test_std <- lapply(xList_test, function(foo) as.data.frame(scale(foo)))

  
  #########################################################################
  ##################### Prediction using posterior samples ################
  #########################################################################
  
  if(bb) {
    message("Prediction using posterior samples...")

    y_samples <- matrix(NA, length(object$errVar_samples), nrow(xList_test[[1]]))
    
    ## Augmented test data
    dataAug_test <- data.Augment(y = NULL, xList_test_std, object$rho_postmed)
    x_aug_test <- as.matrix(dataAug_test$x_aug)
    
    y_aug_samples <- do.call("rbind", lapply(1:nrow(y_samples),
                                             function(tt) {drop(x_aug_test %*% object$beta_samples[tt, ]) + rnorm(nrow(x_aug_test), 0, sqrt(object$errVar_samples[tt]))}))
    y_samples <- y_aug_samples[, 1:nrow(xList_test[[1]])]
    y_pred <- apply(y_samples, 2, median)
    
    ## Returning predictions
    output$y_samples <- y_samples
    output$y_pred <- y_pred + mean(y_test)
    output$y_valid <- y_test
    
  } else {
    message("Prediction using MAP estimator...")
    
    ## Augmented test data
    dataAug_test <- data.Augment(y = NULL, xList_test_std, object$rho_MAP)
    x_aug_test <- as.matrix(dataAug_test$x_aug)
    y_aug_hat <- x_aug_test %*% object$beta_MAP
    y_pred <- y_aug_hat[1:nrow(xList_test[[1]])]
    
    output$y_pred <- y_pred + mean(y_test)
    output$y_valid <- y_test
    
  }
  
  #########################################################################
  ################################ Return #################################
  #########################################################################
  output
}
