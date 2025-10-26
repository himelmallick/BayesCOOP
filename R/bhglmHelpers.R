#' Internal helper functions for BayesCOOP
#'
#' These helpers prepare multimodal data, filter features, construct augmented
#' design matrices for cooperative/consensus-style fusion, and fit weighted
#' spike-and-slab lasso models used inside \code{BayesCOOP()} for Bayesian
#' bootstrap inference and MAP estimation.
#'
#' They are not user-facing and are called internally by \code{BayesCOOP()}.
#'
#' @keywords internal
#' @name bhglmHelpers
#'
#' @importFrom dplyr filter lag
#' @importFrom glmnet glmnet
#' @importFrom MCMCpack rdirichlet
#' @importFrom rmutil rlaplace
#' @importFrom stats rnorm rgamma quantile median
#' @importFrom truncnorm rtruncnorm
#' @importFrom survival cluster
NULL

bhglm.Grouping <- function (all.var, group) {
    n.vars <- length(all.var)
    group.vars <- list()
    if (is.matrix(group)) {
        if (nrow(group) != ncol(group) | ncol(group) > n.vars) 
            stop("wrong dimension for 'group'")
        if (any(rownames(group) != colnames(group))) 
            stop("rownames should be the same as colnames")
        if (any(!colnames(group) %in% all.var)) 
            stop("variabe names in 'group' not in the model predictors")
        group.vars <- colnames(group)
        group <- abs(group)
        wcol <- rowSums(group) - diag(group)
        group <- group/wcol
    }
    else {
        if (is.list(group)) 
            group.vars <- group
        else {
            if (is.numeric(group) & length(group) > 1) {
                group <- sort(group)
                if (group[length(group)] > n.vars) 
                    stop("wrong grouping")
            }
            if (is.numeric(group) & length(group) == 1) 
                group <- as.integer(seq(0, n.vars, length.out = n.vars/group + 
                                            1))
            if (is.null(group)) 
                group <- c(0, n.vars)
            group <- unique(group)
            for (j in 1:(length(group) - 1)) group.vars[[j]] <- all.var[(group[j] + 
                                                                             1):group[j + 1]]
        }
    }
    all.group.vars <- unique(unlist(group.vars))
    if (length(all.group.vars) == n.vars) 
        ungroup.vars <- NULL
    else ungroup.vars <- all.var[which(!all.var %in% all.group.vars)]
    group.new <- c(length(ungroup.vars), length(ungroup.vars) + 
                       cumsum(lapply(group.vars, length)))
    var.new <- c(ungroup.vars, unlist(group.vars))
    list(group = group, group.vars = group.vars, ungroup.vars = ungroup.vars, 
         group.new = group.new, var.new = var.new)
}


bhglm.prepare <- function (x, intercept, prior.mean, prior.sd, prior.scale, prior.df, 
                     group) {
    x0 <- x
    if (intercept) 
        x0 <- x[, -1, drop = FALSE]
    g <- bhglm.Grouping(all.var = colnames(x0), group = group)
    group <- g$group
    group.vars <- g$group.vars
    ungroup.vars <- g$ungroup.vars
    covars <- g$ungroup.vars
    if (is.list(group)) {
        if (length(unlist(group)) > length(unique(unlist(group)))) {
            x1 <- as.data.frame(x0)
            x1 <- x1[, c(covars, unlist(group))]
            g <- c(length(ungroup.vars), length(ungroup.vars) + 
                       cumsum(lapply(group, length)))
            for (j in 1:(length(group) - 1)) group.vars[[j]] <- colnames(x1[, 
                                                                            (g[j] + 1):g[j + 1]])
            x1 <- as.matrix(x1)
            x <- x1
            if (intercept) {
                x <- cbind(1, x)
                colnames(x)[1] <- "(Intercept)"
            }
        }
    }
    J <- NCOL(x)
    if (intercept & J > 1) {
        prior.mean <- c(0, prior.mean)
        prior.scale <- c(prior.scale[1], prior.scale)
        prior.df <- c(prior.df[1], prior.df)
    }
    if (length(prior.mean) < J) 
        prior.mean <- c(prior.mean, rep(prior.mean[length(prior.mean)], 
                                        J - length(prior.mean)))
    if (length(prior.scale) < J) 
        prior.scale <- c(prior.scale, rep(prior.scale[length(prior.scale)], 
                                          J - length(prior.scale)))
    if (length(prior.df) < J) 
        prior.df <- c(prior.df, rep(prior.df[length(prior.df)], 
                                    J - length(prior.df)))
    prior.mean <- prior.mean[1:J]
    prior.scale <- prior.scale[1:J]
    prior.df <- prior.df[1:J]
    prior.df <- ifelse(prior.df == Inf, 1e+10, prior.df)
    if (is.null(prior.sd)) 
        prior.sd <- prior.scale + 0.2
    if (length(prior.sd) < J) 
        prior.sd <- c(prior.sd, rep(prior.sd[length(prior.sd)], 
                                    J - length(prior.sd)))
    prior.sd <- prior.sd[1:J]
    sd.x <- apply(x, 2, sd, na.rm = TRUE)
    min.x.sd <- 1e-04
    prior.sd <- ifelse(sd.x < min.x.sd, 1e-04, prior.sd)
    if (intercept) 
        prior.sd[1] <- 1e+10
    names(prior.mean) <- names(prior.scale) <- names(prior.df) <- names(prior.sd) <- colnames(x)
    if (intercept) 
        covars <- c(colnames(x)[1], covars)
    if (!is.null(covars)) 
        prior.mean[covars] <- 0
    list(x = x, prior.mean = prior.mean, prior.sd = prior.sd, 
         prior.scale = prior.scale, prior.df = prior.df, sd.x = sd.x, 
         min.x.sd = min.x.sd, group = group, group.vars = group.vars, 
         ungroup.vars = ungroup.vars)
}


bhglm.update.ptheta.group <- function (group.vars, p, w, b) {
    f <- function(theta, w, p, bb) {
        sum(p * log(w * theta) + (1 - p) * log(1 - w * theta)) + 
            mean((bb - 1) * log(1 - theta))
    }
    theta <- p
    for (j in 1:length(group.vars)) {
        vars <- group.vars[[j]]
        theta[vars] <- optimize(f, interval = c(0, 1), w = w[vars], 
                                p = p[vars], bb = b[j], maximum = T)$maximum
    }
    theta <- ifelse(theta < 0.01, 0.01, theta)
    theta <- ifelse(theta > 0.99, 0.99, theta)
    theta <- w * theta
    theta
}

bhglm.update.ptheta.network <- function (theta, p, w) {
    phi <- 2
    for (j in 1:length(theta)) {
        mu <- w %*% theta
        m <- mu[j] - w[j, j] * theta[j]
        a <- m * phi
        b <- (1 - m) * phi
        theta[j] <- (p[j] + a)/(1 + a + b)
    }
    theta <- ifelse(theta < 0.01, 0.01, theta)
    theta <- ifelse(theta > 0.99, 0.99, theta)
    theta
}


bhglm.update.scale.p <- function (prior = "mde", df = 1, b0, ss, theta) {
    if (prior == "mde") {
        den0 <- (2 * ss[1])^(-1) * exp(-abs(b0)/ss[1])
        den1 <- (2 * ss[2])^(-1) * exp(-abs(b0)/ss[2])
    }
    if (prior == "mt") {
        den0 <- (ss[1])^(-1) * (1 + b0^2/(df * ss[1]^2))^(-(df + 
                                                                1)/2)
        den1 <- (ss[2])^(-1) * (1 + b0^2/(df * ss[2]^2))^(-(df + 
                                                                1)/2)
    }
    p <- theta * den1/(theta * den1 + (1 - theta) * den0 + 1e-10)
    scale <- 1/((1 - p)/ss[1] + p/ss[2] + 1e-10)
    list(scale = scale, p = p)
}


bhglm.update.theta.weights <- function (eff.hierarchy, eff.parents, p) {
    if (is.null(eff.parents)) 
        stop("'eff.parents' should be given")
    if (!is.list(eff.parents)) 
        stop("'eff.parents' should be a list")
    if (eff.hierarchy == "strong") 
        ww <- lapply(eff.parents, function(x, p) {
            prod(p[x])
        }, p)
    if (eff.hierarchy == "weak") 
        ww <- lapply(eff.parents, function(x, p) {
            mean(p[x])
        }, p)
    ww <- unlist(ww)
    names(ww) <- names(p)
    ww
}

bhglm.covariates <- function (x.con, x.cat, con.rescale = TRUE, cat.center = FALSE, 
                              fill.missing = TRUE, ind.group = NULL) {
    if (missing(x.con)) 
        x.con <- NULL
    if (missing(x.cat)) 
        x.cat <- NULL
    x1 <- x.con
    x2 <- x.cat
    if (!is.null(x1)) {
        x1 <- as.matrix(x1)
        if (is.null(colnames(x1))) {
            warning("some variables have no names.", call. = FALSE)
            colnames(x1) <- paste("x.con", 1:ncol(x1), sep = "")
        }
    }
    if (!is.null(x2)) {
        x2 <- as.matrix(x2)
        if (is.null(colnames(x2))) {
            warning("some variables have no names.", call. = FALSE)
            colnames(x2) <- paste("x.cat", 1:ncol(x2), sep = "")
        }
    }
    x <- cbind(x1, x2)
    func1 <- function(d) {
        na.prop <- length(which(is.na(d)))/length(d)
        return(na.prop)
    }
    na.prop <- apply(x, 2, func1)
    if (any(na.prop > 20/100)) {
        d <- which(na.prop > 20/100)
        warning("the following ", length(d), " variables have more than 20% missing values:", 
                call. = FALSE)
        for (j in 1:length(d)) warning(names(d)[j], call. = FALSE)
    }
    if (!is.null(x1) & con.rescale) 
        x1 <- scale(x1)
    if (!is.null(x2)) {
        X2 <- NULL
        for (i in 1:ncol(x2)) {
            x <- x.new <- x2[, i]
            x <- as.factor(x)
            f <- ~x - 1
            mf <- model.frame(f, na.action = na.pass)
            x.new <- model.matrix(f, mf)[, -1, drop = FALSE]
            colnames(x.new) <- paste(colnames(x2)[i], levels(x)[-1], 
                                     sep = "_")
            if (ncol(x.new) > 1) {
                rsums <- rowSums(x.new)
                if (all(rsums[!is.na(rsums)] == 1)) 
                    x.new <- x.new[, -1, drop = FALSE]
            }
            X2 <- cbind(X2, x.new)
        }
        if (cat.center) {
            x2 <- scale(X2, scale = FALSE)
        }
        else x2 <- X2
    }
    x <- cbind(x1, x2)
    w <- apply(x, 2, var, na.rm = TRUE)
    if (length(is.na(w))) {
        w <- w[!is.na(w)]
        x <- x[, names(w), drop = FALSE]
    }
    if (length(which(w == 0))) {
        x <- x[, which(w != 0), drop = FALSE]
        d <- which(w == 0)
        warning(length(d), " variables with no-variation are removed!", 
                call. = FALSE)
    }
    if (fill.missing & any(is.na(x))) {
        if (!is.null(ind.group)) {
            ind.group <- as.factor(ind.group)
            if (nrow(x) != length(ind.group)) {
                warning("x and ind.group have different obs. Cannot use group information!", 
                        call. = FALSE, immediate. = TRUE)
                ind.group <- NULL
            }
        }
        warning(sum(is.na(x)), " missing values have been filled!", 
                call. = FALSE)
        if (is.null(ind.group)) {
            func2 = function(w) {
                na.index <- which(is.na(w))
                w[na.index] <- mean(w, na.rm = TRUE)
                return(w)
            }
            x <- apply(x, 2, func2)
        }
        if (!is.null(ind.group)) {
            func3 <- function(w) {
                group.means <- tapply(w, ind.group, mean, na.rm = TRUE)
                group.means <- ifelse(is.na(group.means), mean(w, 
                                                               na.rm = TRUE), group.means)
                for (k in 1:length(group.means)) {
                    na.index <- which(is.na(w) & ind.group == names(group.means)[k])
                    w[na.index] <- mean(w[ind.group == names(group.means)[k]], 
                                        na.rm = TRUE)
                }
                return(w)
            }
            x <- apply(x, 2, func3)
        }
    }
    X <- data.frame(x)
    colnames(X) <- colnames(x)
    X
}

bhglm.bmlasso.fit <- function (x, y, family = "gaussian", offset = NULL, epsilon = 1e-04, 
          maxit = 50, init = rep(0, ncol(x)), alpha = 1, ss = c(0.04, 
                                                                0.5), b = 1, group = NULL, theta.weights = NULL, inter.hierarchy = NULL, 
          inter.parents = NULL, Warning = FALSE) 
{
    ss <- sort(ss)
    ss <- ifelse(ss <= 0, 0.001, ss)
    prior.scale <- ss[length(ss)]
    if (family == "cox") 
        intercept <- FALSE
    else intercept <- TRUE
    x0 <- x
    if (intercept) 
        x0 <- cbind(1, x)
    d <- bhglm.prepare(x = x0, intercept = intercept, prior.mean = 0, 
                 prior.sd = 1, prior.scale = prior.scale, prior.df = 1, 
                 group = group)
    x <- d$x
    prior.scale <- d$prior.scale
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
                        alpha = alpha0, lambda = 1/(nrow(x) * ps), standardize = TRUE)
            b <- as.numeric(f$beta)
            if (any(b != 0)) 
                break
        }
    }
    else b <- as.numeric(init)
    names(b) <- colnames(x)
    b <- ifelse(b == 0, 0.001, b)
    init <- b
    devold <- 0
    conv <- FALSE
    for (iter in 1:maxit) {
        if (alpha == 1) 
            out <- bhglm.update.scale.p(prior = "mde", b0 = b[gvars], 
                                  ss = ss, theta = theta)
        else out <- bhglm.update.scale.p(prior = "mt", df = 1e+10, 
                                   b0 = b[gvars], ss = ss, theta = theta)
        prior.scale[gvars] <- out[[1]]
        p <- out[[2]]
        if (!is.matrix(group)) 
            theta <- bhglm.update.ptheta.group(group.vars = group.vars, 
                                         p = p, w = theta.weights, b = bb)
        else theta <- bhglm.update.ptheta.network(theta = theta, p = p, 
                                            w = group)
        if (!is.null(inter.hierarchy)) 
            # theta.weights <- bhglm.update.theta.weights(gvars = gvars, 
            #                                       theta.weights = theta.weights, inter.hierarchy = inter.hierarchy, 
            #                                       inter.parents = inter.parents, p = p)
            theta.weights <- bhglm.update.theta.weights(eff.hierarchy = inter.hierarchy,
                                                        eff.parents = inter.parents, p = p)
        Pf <- 1/(prior.scale + 1e-10)
        f <- glmnet(x = x, y = y, family = family, offset = offset, 
                    alpha = alpha, penalty.factor = Pf, lambda = sum(Pf)/(nrow(x) * 
                                                                              ncol(x)), standardize = FALSE)
        b <- as.numeric(f$beta)
        names(b) <- colnames(x)
        dev <- deviance(f)
        if (abs(dev - devold)/(0.1 + abs(dev)) < epsilon & iter > 
            5) {
            conv <- TRUE
            break
        }
        else devold <- dev
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
    #     f$dispersion <- bglm(y ~ f$linear.predictors - 1, start = 1, 
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


bhglm.bmlasso <- function (x, y, family = c("gaussian", "binomial", "poisson", 
                           "cox"), offset = NULL, epsilon = 1e-04, maxit = 50, init = NULL, 
          alpha = c(1, 0), ss = c(0.04, 0.5), b = 1, group = NULL, 
          theta.weights = NULL, inter.hierarchy = NULL, inter.parents = NULL, 
          Warning = FALSE, verbose = FALSE) 
{
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
    f <- bhglm.bmlasso.fit(x = x, y = y, family = family, offset = offset, 
                     epsilon = epsilon, maxit = maxit, init = init, group = group, 
                     alpha = alpha, ss = ss, b = b, theta.weights = theta.weights, 
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



