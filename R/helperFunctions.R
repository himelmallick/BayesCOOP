library(BhGLM)
library(caret)
library(tidyverse)

###############################################
############## helper functions ###############
###############################################

# Function to scale a matrix
scalex <- function(x, cmeans, csds) {
    ## x is a matrix
    x2 <- matrix(NA, nrow(x), ncol(x))
    for(i in 1:ncol(x)) {
        x2[, i] <- (x[, i] - cmeans[i]) / csds[i]
    }
    x2
}


# Function to reformat the data for cooperative learning application
gen_datalist <- function(data) {
    feature_table <- data$feature_table
    sample_metadata <- data$sample_metadata
    feature_metadata <- data$feature_metadata
    
    ## Separate omics layers
    feature_metadata$featureType <- as.factor(feature_metadata$featureType)
    name_layers <- with(droplevels(feature_metadata),
                       list(levels = levels(featureType)),
                       nlevels = nlevels(featureType))$levels
    
    ## Extract y
    y <- as.matrix(sample_metadata$Y)
    
    ## Construct xList
    xList <- vector("list", length = length(name_layers))
    names(xList) <- name_layers
    
    ## Store data matrices
    for (i in seq_along(name_layers)) {
        ## Filter genes which belong to featureType[i]
        include_list <- feature_metadata %>% dplyr::filter(featureType == name_layers[i])
        
        ## Subset feature rows for features in Type[i]
        t_dat_slice <- feature_table[rownames(feature_table) %in% include_list$featureID, ]
        
        ## Assign feature rows for Type[i] to dataList[i]
        xList[[i]] <- as.data.frame(t(t_dat_slice))
        rm(t_dat_slice); rm(include_list)
    }
    
    list(xList = xList, y = y)
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
    y_aug <- c(y, rep(0, length(y)))
    tdataList <- list()
    tdataList[[1]] <- rbind.data.frame(dataList[[1]], -sqrt(rho) * dataList[[1]])
    tdataList[[2]] <- rbind.data.frame(dataList[[2]], sqrt(rho) * dataList[[2]])
    x_aug <- do.call("cbind", tdataList)
    return(list(y_aug = y_aug, x_aug = x_aug))
}


data.Augment.mviews <- function(y, dataList, rho){
    p_x <- lapply(dataList, ncol)
    nviews <- length(dataList)
    pairs <- apply(utils::combn(nviews, 2), 2, identity, simplify = FALSE)
    npairs <- length(pairs)
    x <- do.call(cbind, dataList)
    
    rows <- lapply(pairs, make_row, x_list = dataList, p_x = p_x, rho = rho)
    # matching the column names
    for(i in 1:length(rows)){
        colnames(rows[[i]]) <- colnames(x)
    }
    xt <- do.call(rbind, c(list(x), rows))
    yt <- c(y, rep(0, length(pairs) * length(y)))
    return(list(y_aug = yt, x_aug = xt))
}


# Function to filter features according to a fixed abundance and prevalence threshold
filter_features <- function(x, abd_threshold = 0, prev_threshold = 0.1) {
    # x is a data frame
    x <- x[, colMeans(x > abd_threshold) > prev_threshold]
    nzv_x <- nearZeroVar(x, names = TRUE)
    x <- x[, setdiff(names(x), nzv_x)] 
    return(x)
}

## bmlasso.fit with likelihood weights for Bayesian Bootstrap approach

# jitter: a location-shift vector to create a random shift in the centering of the prior of beta. {cite Nie & Rockova, 2022, JASA paper}.
# lwood_weights: randomly samples weights for re-weighting the likelihood function.
# Check BhGLM package for other arguments
bmlasso.fit.weighted <- function (x, y, family = "gaussian", offset = NULL, epsilon = 1e-04, 
                                  maxit = 50, init = rep(0, ncol(x)), alpha = 1,
                                  ss = c(0.05, 1), b = 1, lhood_weights = lhood_weights, 
                                  jitter = 0, group = NULL, theta.weights = NULL, inter.hierarchy = NULL, 
                                  inter.parents = NULL, Warning = FALSE) 
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
    
    d <- BhGLM:::prepare(x = x0, intercept = intercept, prior.mean = jitter,
                 prior.sd = 1, prior.scale = prior.scale, prior.df = 1,
                 group = group)
    x <- d$x
    prior.scale <- d$prior.scale
    
    ###############################################
    prior.mean = d$prior.mean[-length(d$prior.mean)] # Additional line
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
    b <- b[1:length(group.vars)]
    bb <- b
    if (is.null(init)) {
        for (k in 1:5) {
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
    names(b) = gvars                               # names(b) <- colnames(x)
    b <- ifelse(b == 0, 0.001, b)
    init <- b
    devold <- 0
    conv <- FALSE
    for (iter in 1:maxit) {
        if (alpha == 1) 
            out <- BhGLM:::update.scale.p(prior = "mde", b0 = (b[gvars] - prior.mean[-1]), 
                                  ss = ss, theta = theta) else 
                                      out <- BhGLM:::update.scale.p(prior = "mt", df = 1e+10, 
                                   b0 = b[gvars], ss = ss, theta = theta)
        # prior.scale[gvars] <- out[[1]]
        prior.scale = out[[1]]
        p <- out[[2]]
        if (!is.matrix(group)) 
            theta <- BhGLM:::update.ptheta.group(group.vars = group.vars, 
                                         p = p, w = theta.weights, b = bb) else 
                                             theta <- BhGLM:::update.ptheta.network(theta = theta, p = p, 
                                            w = group)
        if (!is.null(inter.hierarchy)) 
            theta.weights <- BhGLM:::update.theta.weights(gvars = gvars, 
                                                  theta.weights = theta.weights, inter.hierarchy = inter.hierarchy, 
                                                  inter.parents = inter.parents, p = p)
        Pf <- 1/(prior.scale + 1e-10)
        f <- glmnet(x = x, y = y, family = family, offset = offset, 
                    alpha = alpha, penalty.factor = Pf, weights = lhood_weights,
                    lambda = sum(Pf)/(nrow(x) * ncol(x)), standardize = FALSE)
        b <- as.numeric(f$beta)
        names(b) = gvars                   # names(b) <- colnames(x)
        dev <- deviance(f)
        if (abs(dev - devold)/(0.1 + abs(dev)) < epsilon & iter > 
            5) {
            conv <- TRUE
            break
        } else devold <- dev
    }
    if (Warning & !conv) 
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
    #     f$dispersion <- BhGLM:::bglm(y ~ f$linear.predictors - 1, start = 1, 
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

bmlasso.weighted <- function (x, y, family = c("gaussian", "binomial", "poisson", 
                           "cox"), offset = NULL, epsilon = 1e-04, maxit = 50, init = NULL, 
          alpha = c(1, 0), ss = c(0.04, 0.5), b = 1, group = NULL, lhood_weights = lhood_weights, jitter = NULL,
          theta.weights = NULL, inter.hierarchy = NULL, inter.parents = NULL, 
          Warning = FALSE, verbose = FALSE) 
{
    if (!requireNamespace("glmnet")) 
        install.packages("glmnet")
    require(glmnet)
    start.time <- Sys.time()
    call <- match.call()
    x <- as.matrix(x)
    if (is.null(colnames(x))) 
        colnames(x) <- paste("x", 1:ncol(x), sep = "")
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
                     Warning = Warning)
    f$call <- call
    if (family == "cox") 
        class(f) <- c(class(f), "bmlasso", "COXPH")
    else class(f) <- c(class(f), "bmlasso", "GLM")
    stop.time <- Sys.time()
    minutes <- round(difftime(stop.time, start.time, units = "min"), 
                     3)
    if (verbose) {
        cat("EM Coordinate Decent Iterations:", f$iter, "\n")
        cat("Computational time:", minutes, "minutes \n")
    }
    return(f)
}


