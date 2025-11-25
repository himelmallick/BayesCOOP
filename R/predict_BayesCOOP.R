#' @title predict_BayesCOOP
#'
#' @description This function predicts the response from a fitted object using BayesCOOP methodology. Currently, only continuous outcomes are supported.
#'
#' @param object the object returned by the function \code{\link{BayesCOOP}} implemented on the training dataset.
#' @param newdata a list of feature_table, sample_metadata and feature_metadata from validation dataset (unstandardized). See \emph{IntegratedLearner} for more details.
#' @param family currently supports only Gaussian family. Default: "gaussian".
#' @param bb logical. If TRUE, run full Bayesian bootstrap inference; if FALSE, run MAP
#' estimation over a user-supplied rho grid in \code{control}. Default: TRUE.
#' @param filter logical. If TRUE, apply abundance/prevalence filtering to features. Default: TRUE.
#' @param abd_thresh minimum abundance threshold for keeping a feature. Default: 0.
#' @param prev_thresh minimum prevalence threshold (proportion of samples above abd_thresh).
#' Default: 0.1.
#' @param Warning logical. If TRUE, emit non-convergence warnings. Default: TRUE.
#' @param verbose logical. If TRUE, print iteration counts and runtime. Default: TRUE.
#' @param control a named list with element \code{rho}, giving one or more candidate rho values
#' to try when \code{bb = FALSE}.
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
#'
#' @export
predict_BayesCOOP <- function(object, newdata, family = "gaussian",
                      bb = TRUE, filter = TRUE, abd_thresh = 0, prev_thresh = 0.1,
                      Warning = TRUE, verbose = TRUE, control = list()) {

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

    ## Filtering features
    if(filter == TRUE){
        xList_test <- lapply(xList_test, function(foo) filter_features(foo, abd_thresh, prev_thresh))
    }

    ## Considering the common features between train and test set
    keep_features <- vector("list", length = length(object$feature_names))
    for(i in 1:length(keep_features)) {
        keep_features[[i]] <- intersect(object$feature_names[[i]], names(xList_test[[i]]))
        xList_test[[i]] <- as.matrix(xList_test[[i]][, keep_features[[i]], drop = FALSE])
    }

    ## Scaling the predictors
    xList_test_std <- lapply(xList_test, function(foo) as.data.frame(scale(foo)))

    ## centering the response
    #y_test_ctd <- y_test - mean(y_test)

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
        output$y_pred <- y_pred
        output$y_valid <- y_test

    } else {
        message("Prediction using MAP estimator...")

        ## Augmented test data
        dataAug_test <- data.Augment(y = NULL, xList_test_std, object$rho_MAP)
        x_aug_test <- as.matrix(dataAug_test$x_aug)
        y_aug_hat <- x_aug_test %*% object$beta_MAP
        y_pred <- y_aug_hat[1:nrow(xList_test[[1]])]

        output$y_pred <- y_pred
        output$y_valid <- y_test

    }

    #########################################################################
    ################################ Return #################################
    #########################################################################
    output
}
