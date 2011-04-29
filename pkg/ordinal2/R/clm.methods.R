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
  
logLik.clm <- function(object, ...)
  structure(object$logLik, df = object$edf, class = "logLik")

extractAIC.clm <- function(fit, scale = 0, k = 2, ...) {
  edf <- fit$edf
  c(edf, -2*fit$logLik + k * edf)
}

### NOTE: AIC.clm implicitly defined via logLik.clm

anova.clm <- function (object, ..., test = c("Chisq", "none"))
### make more appropriate printing of model formulae, potential
### threshold structures and link functions.
{
  test <- match.arg(test)
  dots <- list(...)
  if (length(dots) == 0)
    stop('anova is not implemented for a single "clm" object')
  mlist <- list(object, ...)
  nt <- length(mlist)
  dflis <- sapply(mlist, function(x) x$df.residual)
  s <- order(dflis, decreasing = TRUE)
  mlist <- mlist[s]
  if (any(!sapply(mlist, inherits, "clm")))
    stop('not all objects are of class "clm"')
  ns <- sapply(mlist, function(x) length(x$fitted.values))
  if(any(ns != ns[1]))
    stop("models were not all fitted to the same dataset")
  rsp <- unique(sapply(mlist, function(x) {
                       tmp <- attr(x$model, "terms")
                       class(tmp) <- "formula"
                       paste(tmp[2]) } ))
  mds <- sapply(mlist, function(x) {
      tmp1 <- attr(x$model, "terms")
      class(tmp1) <- "formula"
      paste(tmp1[3]) })
### FIXME: Extract formulae from the matched call instead to get the
### random effects included as well?
  dfs <- dflis[s]
  lls <- sapply(mlist, function(x) -2*x$logLik)
  tss <- c("", paste(1:(nt - 1), 2:nt, sep = " vs "))
  df <- c(NA, -diff(dfs))
  x2 <- c(NA, -diff(lls))
  pr <- c(NA, 1 - pchisq(x2[-1], df[-1]))
  out <- data.frame(Model = mds, Resid.df = dfs, '-2logLik' = lls,
                    Test = tss, Df = df, LRtest = x2, Prob = pr)
  names(out) <- c("Model", "Resid. df", "-2logLik", "Test",
                  "   Df", "LR stat.", "Pr(Chi)")
  if (test == "none") out <- out[, 1:6]
  class(out) <- c("Anova", "data.frame")
  attr(out, "heading") <-
    c("Likelihood ratio tests of cumulative link models\n",
      paste("Response:", rsp))
  out
}

