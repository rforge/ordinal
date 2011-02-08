print.clmm <- function(x, ...)
{
  cat("Cumulative Link Mixed Model fitted with the Laplace approximation\n",
      fill=TRUE)
  if(!is.null(cl <- x$call)) {
    cat("Call:\n")
    dput(cl, control=NULL)
  }

  cat("\nRandom effects:\n")
  varMat <- matrix(c(x$stDev^2, x$stDev), nrow =
                   length(x$stDev), ncol=2)
  rownames(varMat) <- names(x$stDev)
  colnames(varMat) <- c("Var", "Std.Dev")
  print(varMat, ...)
  
  if(length(x$beta)) {
    cat("\nCoefficients:\n")
    print(x$beta, ...)
  } else 
    cat("\nNo Coefficients\n")
  
  if(length(x$alpha) > 0) {
    cat("\nThresholds:\n")
    print(x$alpha, ...)
  }
  cat("\nlog-likelihood:", format(x$logLik, nsmall=2), "\n")
  cat("AIC:", format(-2*x$logLik + 2*x$edf, nsmall=2), "\n")
  if(nzchar(mess <- naprint(x$na.action))) cat("(", mess, ")\n", sep="")
  invisible(x)
}

vcov.clmm <- function(object, ...)
{
    if(is.null(object$Hessian)) {
        stop("Model needs to be fitted with Hess = TRUE")
    }
    dn <- dimnames(object$Hessian)
    structure(ginv(object$Hessian), dimnames = dn)
}

summary.clmm <- function(object, digits = max(3, .Options$digits - 3),
                         correlation = FALSE, ...)
{
  object$varMat <- matrix(c(object$stDev^2, object$stDev),
                          nrow = length(object$stDev), ncol=2)
  rownames(object$varMat) <- names(object$stDev)
  colnames(object$varMat) <- c("Var", "Std.Dev")
  
  npar <- length(object$alpha) + length(object$beta)
  coef <- matrix(0, npar, 4,
                 dimnames = list(names(object$coefficients[1:npar]),
                   c("Estimate", "Std. Error", "z value", "Pr(>|z|)")))
  coef[, 1] <- object$coefficients[1:npar]
  vc <- try(vcov(object), silent = TRUE)
  if(class(vc) == "try-error") {
    warning("Variance-covariance matrix of the parameters is not defined")
    coef[, 2:4] <- NaN
    if(correlation) warning("Correlation matrix is unavailable")
    object$condHess <- NaN
  }
  else {
    coef[, 2] <- sd <- sqrt(diag(vc)[1:npar])
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
  class(object) <- "summary.clmm"
  object
}

print.summary.clmm <- function(x, digits = x$digits, signif.stars =
                              getOption("show.signif.stars"), ...)
{
### FIXME: what about the the signif.stars ???
### FIXME: print lmer-like information about data size, logLik, AIC,
### BIC, etc.
  cat("Cumulative Link Mixed Model fitted with the Laplace approximation\n",
      fill=TRUE)
  if(!is.null(cl <- x$call)) {
    cat("Call:\n")
    dput(cl, control=NULL)
  }
  
  cat("\nRandom effects:\n")
  print(x$varMat, ...)
  
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
