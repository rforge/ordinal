print.clm <- function(x, ...)
{
  cat("formula:", deparse(x$call$formula), fill=TRUE)
  if(!is.null(data.name <- x$call$data))
    cat("data:   ", deparse(data.name), fill=TRUE)
  if(!is.null(x$call$subset))
    cat("subset: ", deparse(x$call$subset), fill=TRUE)
  cat("\n")

  print(x$info, row.names=FALSE, right=FALSE)

  if(length(x$beta)) {
    cat("\nCoefficients:\n")
    print(x$beta, ...)
  } else {
    cat("\nNo Coefficients\n")
  }
  if(length(x$alpha) > 0) {
    cat("\nThresholds:\n")
    print(x$alpha, ...)
  }

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
  if(is.null(object$Hessian))
    stop("Model needs to be fitted with Hess = TRUE")
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
  object$info$cond.H <- formatC(object$condHess, digits=1, format="e")
  object$coefficients <- coef
  object$digits <- digits
  class(object) <- "summary.clm"
  object
}

print.summary.clm <- function(x, digits = x$digits, signif.stars =
                              getOption("show.signif.stars"), ...)
{
  cat("formula:", deparse(x$call$formula), fill=TRUE)
  if(!is.null(data.name <- x$call$data))
    cat("data:   ", deparse(data.name), fill=TRUE)
  if(!is.null(x$call$subset))
    cat("subset: ", deparse(x$call$subset), fill=TRUE)
  cat("\n")

  print(x$info, row.names=FALSE, right=FALSE)
  
  nbeta <- length(x$beta)
  nalpha <- length(x$alpha)
  if(nbeta > 0) {
    cat("\nCoefficients:\n")
    printCoefmat(x$coefficients[nalpha + 1:nbeta, , drop=FALSE],
                 digits=digits, signif.stars=signif.stars,
                 has.Pvalue=TRUE, ...) 
  } else {
    cat("\nNo Coefficients\n")
  }
  if(nalpha > 0) { ## always true
    cat("\nThreshold coefficients:\n")
    printCoefmat(x$coefficients[seq_len(nalpha), -4, drop=FALSE],
                 digits=digits, has.Pvalue=FALSE, signif.stars=FALSE,
                 ...) 
  }

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

## anova.clm <- function (object, ..., test = c("Chisq", "none"))
## ### make more appropriate printing of model formulae, potential
## ### threshold structures and link functions.
## {
##   test <- match.arg(test)
##   dots <- list(...)
##   if (length(dots) == 0)
##     stop('anova is not implemented for a single "clm" object')
##   mlist <- list(object, ...)
##   nt <- length(mlist)
##   dflis <- sapply(mlist, function(x) x$df.residual)
##   s <- order(dflis, decreasing = TRUE)
##   mlist <- mlist[s]
##   if (any(!sapply(mlist, inherits, "clm")))
##     stop('not all objects are of class "clm"')
##   ns <- sapply(mlist, function(x) length(x$fitted.values))
##   if(any(ns != ns[1]))
##     stop("models were not all fitted to the same dataset")
##   rsp <- unique(sapply(mlist, function(x) {
##                        tmp <- attr(x$model, "terms")
##                        class(tmp) <- "formula"
##                        paste(tmp[2]) } ))
##   mds <- sapply(mlist, function(x) {
##       tmp1 <- attr(x$model, "terms")
##       class(tmp1) <- "formula"
##       paste(tmp1[3]) })
## ### FIXME: Extract formulae from the matched call instead to get the
## ### random effects included as well?
##   dfs <- dflis[s]
##   lls <- sapply(mlist, function(x) -2*x$logLik)
##   tss <- c("", paste(1:(nt - 1), 2:nt, sep = " vs "))
##   df <- c(NA, -diff(dfs))
##   x2 <- c(NA, -diff(lls))
##   pr <- c(NA, 1 - pchisq(x2[-1], df[-1]))
##   out <- data.frame(Model = mds, Resid.df = dfs, '-2logLik' = lls,
##                     Test = tss, Df = df, LRtest = x2, Prob = pr)
##   names(out) <- c("Model", "Resid. df", "-2logLik", "Test",
##                   "   Df", "LR stat.", "Pr(Chi)")
##   if (test == "none") out <- out[, 1:6]
##   class(out) <- c("Anova", "data.frame")
##   attr(out, "heading") <-
##     c("Likelihood ratio tests of cumulative link models\n",
##       paste("Response:", rsp))
##   out
## }

anova.clm <- function(object, ...)
### requires that clm objects have components:
###  edf: no. parameters used
###  call$formula
###  link (character)
###  threshold (character)
###  logLik
###  
{
  dots <- list(...)
  if (length(dots) == 0)
    stop('anova is not implemented for a single "clm" object')
  mlist <- list(object, ...)
  no.par <- sapply(mlist, function(x) x$edf)
  ## order list with increasing no. par:
  ord <- order(no.par, decreasing=TRUE)
  mlist <- mlist[ord]
  no.par <- no.par[ord]
  no.tests <- length(mlist)
  forms <- sapply(mlist, function(x) deparse(x$call$formula))
  links <- sapply(mlist, function(x) x$link)
  thres <- sapply(mlist, function(x) x$threshold)
  models <- data.frame(forms, links, thres)
  models.names <- c('formula:', "link:", "threshold:")
  AIC <- sapply(mlist, function(x) -2*x$logLik + 2*x$edf)
  logLiks <- sapply(mlist, function(x) x$logLik)
  statistic <- c(NA, -2*diff(sapply(mlist, function(x) x$logLik)))
  df <- c(NA, -diff(no.par))
  pval <- c(NA, pchisq(statistic[-1], df[-1], lower.tail=FALSE))
  pval[!is.na(df) & df==0] <- NA 
  tab <- data.frame(no.par, AIC, logLiks, statistic, df, pval) 
  tab.names <- c("no.par", "AIC", "logLik", "LR.stat", "df",
                 "Pr(>Chisq)")
  colnames(tab) <- tab.names
  rownames(tab) <- rownames(models) <- 1:no.tests
  colnames(models) <- models.names
  attr(tab, "models") <- models
  attr(tab, "heading") <-
    "Likelihood ratio tests of cumulative link models:\n"
  class(tab) <- c("anova.clm", "data.frame")
  tab
}

print.anova.clm <-
  function(x, digits=max(getOption("digits") - 2, 3),
           signif.stars=getOption("show.signif.stars"), ...) 
{
  if (!is.null(heading <- attr(x, "heading"))) 
    cat(heading, "\n")
  models <- attr(x, "models")
  print(models, right=FALSE)
  printCoefmat(x, digits=digits, signif.stars=signif.stars,
               tst.ind=4, cs.ind=NULL, # zap.ind=2, #c(1,5),
               P.values=TRUE, has.Pvalue=TRUE, na.print="")
}

