print.clm <- function(x, ...)
{
  if(!is.null(cl <- x$call)) {
    cat("Call:\n")
    dput(cl, control=NULL)
  }
  if(length(x$beta)) {
    cat("\nCoefficients:\n")
    print(x$beta, ...)
  } else {
    cat("\nNo Coefficients\n")
  }
  if(length(x$alpha) > 0) {
    cat("\nThresholds:\n")
    print(x$alpha, ...)
    ## if(x$threshold != "flexible") {
    ##   cat("\nThresholds:\n")
    ##   print(x$Theta, ...)
    ## }
  }
    cat("\nlog-likelihood:", format(x$logLik, nsmall=2), "\n")
    cat("AIC:", format(-2*x$logLik + 2*x$edf, nsmall=2), "\n")
    if(nzchar(mess <- naprint(x$na.action))) cat("(", mess, ")\n", sep="")
    invisible(x)
}

vcov.clm <- function(object, ...)
{
  dn <- names(object$coefficients)
  H <- object$Hessian
  ## To handle NaNs in the Hessian resulting from parameter
  ## unidentifiability:  
  if(any(His.na <- !is.finite(H))) {
    H[His.na] <- 0
    VCOV <- MASS::ginv(H)
    VCOV[His.na] <- NaN
  }
  else
    VCOV <- MASS::ginv(H)
  structure(VCOV, dimnames = list(dn, dn))
}

summary.clm <- function(object, digits = max(3, .Options$digits - 3),
                        correlation = FALSE, ...)
{
  coef <- matrix(0, object$edf, 4,
                 dimnames = list(names(object$coefficients),
                   c("Estimate", "Std. Error", "z value", "Pr(>|z|)")))
  coef[, 1] <- object$coefficients
  vc <- try(vcov(object), silent = TRUE)
  if(class(vc) == "try-error") {
    warning("Variance-covariance matrix of the parameters is not defined")
    coef[, 2:4] <- NaN
    if(correlation) warning("Correlation matrix is unavailable")
    object$condHess <- NaN
  }
  else {
    coef[, 2] <- sd <- sqrt(diag(vc))
    ## Cond is Inf if Hessian contains NaNs:
    object$condHess <-
      if(any(is.na(object$Hessian))) Inf
      else with(eigen(object$Hessian, only.values = TRUE),
                abs(max(values) / min(values)))
    coef[, 3] <- coef[, 1]/coef[, 2]
    coef[, 4] <- 2 * pnorm(abs(coef[, 3]), lower.tail=FALSE)
    if(correlation)
        object$correlation <-
          (vc / sd) / rep(sd, rep(object$edf, object$edf))
  }
  object$coefficients <- coef
  object$digits <- digits
  class(object) <- "summary.clm"
  object
}

print.summary.clm <- function(x, digits = x$digits, signif.stars =
                              getOption("show.signif.stars"), ...)
{
### FIXME: what about the the signif.stars ???
### FIXME: print lmer-like information about data size, logLik, AIC,
### BIC, etc.
  if(!is.null(cl <- x$call)) {
    cat("Call:\n")
    dput(cl, control=NULL)
  }
  coef <- format(round(x$coefficients, digits=digits))
  coef[,4] <- format.pval(x$coefficients[, 4])
  nbeta <- length(x$beta); nalpha <- length(x$alpha)
  if(nbeta > 0) {
    cat("\nCoefficients:\n")
    print(coef[nalpha + 1:nbeta, , drop=FALSE], quote = FALSE, ...)
  } else {
    cat("\nNo Coefficients\n")
  }
  if(nalpha > 0) { ## always true
    cat("\nThreshold coefficients:\n")
    print(coef[seq_len(nalpha), -4, drop=FALSE], quote = FALSE, ...)
  }
  
  cat("\nlog-likelihood:", format(x$logLik, nsmall=2), "\n")
  cat("AIC:", format(-2*x$logLik + 2*x$edf, nsmall=2), "\n")
  cat("Condition number of Hessian:", format(x$condHess, nsmall=2), "\n")
  if(nzchar(mess <- naprint(x$na.action))) cat("(", mess, ")\n", sep="")
  if(!is.null(correl <- x$correlation)) {
    cat("\nCorrelation of Coefficients:\n")
    ll <- lower.tri(correl)
    correl[ll] <- format(round(correl[ll], digits))
    correl[!ll] <- ""
    print(correl[-1, -ncol(correl)], quote = FALSE, ...)
  }
  invisible(x)
}

slice <- function(object, ...) {
  UseMethod("slice")
}


slice.clm <-
  function(x, which = seq_along(par), lambda = 3, grid = 1e2,
           quad = TRUE)
{
  ## argument matching and testing:
  stopifnot(is.numeric(lambda) && lambda > 0)
  stopifnot(is.numeric(grid) && grid >= 1)
  grid <- as.integer(grid)
  par <- coef(x)
  stopifnot(is.numeric(which))
  stopifnot(length(which) == length(unique(which)))
  stopifnot(all(which >= 1) && all(which <= length(par)))
  which <- as.integer(which)
  parNames <- names(par)
  npar <- length(par)
  ml <- x$logLik
  
  ## get environment corresponding to x:
  rho <- update(x, doFit = FALSE)
  names(par) <- NULL
  rho$par <- par ## set rho$par to mle
  stopifnot(isTRUE(all.equal(clm.nll(rho), -x$logLik)))

  curv <- 1/diag(x$Hess) ## curvature in nll wrt. par
  par.range <- par + sqrt(curv) %o% c(-lambda, lambda)
  ## par.seq - list of length npar:
  par.seq <- sapply(which, function(ind) {
    seq(par.range[ind, 1], par.range[ind, 2], length = grid) }, 
                    simplify = FALSE)
  ## compute relative nll for all par.seq for each par:
  nll <- lapply(which, function(par.ind) { # for each par
    rho$par <- par ## reset par values to MLE
    sapply(par.seq[[par.ind]], function(par.val) { # for each val
      rho$par[par.ind] <- par.val
      clm.nll(rho) + ml ## relative nll
    })
  })
  
  ## collect results in a list of data.frames:
  res <- lapply(which, function(i) {
    structure(data.frame(par.seq[[i]], nll[[i]]),
              names = c(parNames[i], "nll"))
  })

  ## set attributes:
  names(res) <- parNames[which]
  attr(res, "original.fit") <- x
  class(res) <- "slice.clm"

  if(!quad) return(res)
  ## compute quadratic approx to nll:
  Quad <- function(par, mle, curv)
    ((mle - par)^2 / curv / 2)
  for(i in which) 
    res[[i]]$quad <- Quad(par.seq[[i]], par[i], curv[i])
  
  return(res)
}

plot.slice.clm <-
  function(x, which = seq_along(x), type = c("quadratic", "linear"),
           plot.mle = TRUE, ...)
{
  type <- match.arg(type)
  stopifnot(is.numeric(which))
  which <- as.integer(which)
  ml <- attr(x, "original.fit")$logLik

  ## take the signed sqrt of nll and quad:
  if(type == "linear") {
    sgn.sqrt <- function(par, mle, nll)
      (2 * (par > mle) - 1) * sqrt(nll)
    mle <- coef(attr(x, "original.fit"))
    for(i in which) {
      x[[i]]$nll <- sgn.sqrt(x[[i]][1], mle[i], x[[i]]$nll)
      if(!is.null(x[[i]]$quad))
        x[[i]]$quad <- sgn.sqrt(x[[i]][1], mle[i], x[[i]]$quad)
    }
  }
 
  ## actual plotting:
  for(i in which) {
    z <- x[[i]]
    plot(z[1:2], type = "l", ...)
    if(!is.null(z$quad))
      lines(z[[1]], z[[3]], lty = 2)
    if(plot.mle && type == "quadratic")
      abline(v = coef(attr(x, "original.fit"))[i])
  }
  
  return(invisible())
}
  
