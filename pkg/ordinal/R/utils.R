getFitted <- function(eta1, eta2, pfun, ...) {
  ## eta1, eta2: linear predictors
  ## pfun: cumulative distribution function
  ##
  ## Compute fitted values while maintaining high precision in the
  ## result - if eta1 and eta2 are both large, fitted is the
  ## difference between two numbers very close to 1, which leads to
  ## imprecision and potentially errors.
  ##
  ## Note that (eta1 > eta2) always holds, hence (eta2 > 0) happens
  ## relatively rarely.
  k2 <- eta2 > 0
  fitted <- pfun(eta1) - pfun(eta2)
  fitted[k2] <- pfun(eta2[k2], lower.tail=FALSE) -
    pfun(eta1[k2], lower.tail=FALSE)
  fitted
}

getFittedC <-
  function(eta1, eta2,
           link = c("logit", "probit", "cloglog", "loglog", "cauchit",
             "Aranda-Ordaz", "log-gamma"), lambda=1)
### Same as getFitted only this is implemented in C and handles all
### link functions including the flexible ones.
{
  link <- match.arg(link)
  .Call("get_fitted", eta1, eta2, link, lambda)
}

getWeights <- function(mf) {
### mf - model.frame
    n <- nrow(mf)
    if(is.null(wts <- model.weights(mf))) wts <- rep(1, n)
    ## if (any(wts <= 0))
    ##   stop(gettextf("non-positive weights are not allowed"),
    ##        call.=FALSE)
### NOTE: We do not remove observations where weights == 0, because
### that could be a somewhat surprising behaviour. It would also
### require that the model.frame be evaluated all over again to get
### the right response vector with the right number of levels.
    if(length(wts) && length(wts) != n)
        stop(gettextf("number of weights is %d should equal %d (number of observations)",
                      length(wts), n), call.=FALSE)
    if(any(wts < 0))
        stop(gettextf("negative weights are not allowed"),
             call.=FALSE)
    ## if(any(wts == 0)) {
    ##     y <- model.response(mf, "any")
    ##     if(any(table(y[wts > 0]) == 0))
    ##         stop(gettextf("zero positive weights for one or more response categories"),
    ##              call.=FALSE)
    ## }
    return(as.double(wts))
}

getOffset <- function(mf) {
### mf - model.frame
  n <- nrow(mf)
  if(is.null(off <- model.offset(mf))) off <- rep(0, n)
  if(length(off) && length(off) != n)
    stop(gettextf("number of offsets is %d should equal %d (number of observations)",
                  length(off), n), call.=FALSE)
  return(as.double(off))
}

getFullForm <- function(form, ..., envir=parent.frame()) {
### collect terms in several formulas in a single formula
### sets the environment of the resulting formula to envir.
  forms <- list(...)
  if(lf <- length(forms)) {
    rhs <- character(0)
    ## Collect rhs terms in a single vector of rh-sides:
    for(i in 1:lf) {
      rhs <- c(rhs, Deparse(forms[[i]][[2]]))
      if(length(forms[[i]]) >= 3)
        rhs <- c(rhs, Deparse(forms[[i]][[3]]))
    }
    ## add '+' inbetween terms:
    rhs <- paste(rhs, collapse=" + ")
    ## combine if 'deparse(form)' is a (long) vector:
    form2 <- paste(deparse(form, width.cutoff=500L), collapse=" ")
    ## combine form2 and rhs into a single string:
    form <- paste(form2, rhs, sep=" + ")
  }
  return(as.formula(form, env=envir))
}

## getFullForm <- function(form, ..., envir=parent.frame()) {
## ### collect terms in several formulas in a single formula (on the rhs)
## ### sets the environment of the resulting formula to envir.
##   forms <- list(form, ...)
##   allVars <- unlist(sapply(forms, all.vars))
##   rhs <- paste(allVars, collapse=" + ")
##   form <- paste("~", rhs)
##   return(as.formula(form, env=envir))
## }

Trace <- function(iter, stepFactor, val, maxGrad, par, first=FALSE) {
    t1 <- sprintf(" %3d:  %-5e:    %.3f:   %1.3e:  ",
                  iter, stepFactor, val, maxGrad)
    t2 <- formatC(par)
    if(first)
        cat("iter:  step factor:     Value:     max|grad|:   Parameters:\n")
    cat(t1, t2, "\n")
}

response.name <- function(terms) {
  vars <- as.character(attr(terms, "variables"))
  vars[1 + attr(terms, "response")]
}

getB <- function(y, NOM=NULL, X=NULL, offset=NULL, tJac=NULL) {
### FIXME: Is this function ever used?
### NOTE: no tests that arguments conform.
  nlev <- nlevels(y)
  n <- length(y)
  B2 <- 1 * (col(matrix(0, n, nlev)) == c(unclass(y)))
  o1 <- c(1e5 * B2[, nlev]) - offset
  o2 <- c(-1e5 * B2[,1]) - offset
  B1 <- B2[, -(nlev), drop = FALSE]
  B2 <- B2[, -1, drop = FALSE]
  ## adjust B1 and B2 for structured thresholds:
  if(!is.null(tJac)) {
    B1 <- B1 %*% tJac
    B2 <- B2 %*% tJac
  }
  ## update B1 and B2 with nominal effects:
  if(NCOL(NOM) > 1) { ## !is.null(NOM) && ncol(NOM) > 1) {
    ## if !is.null(NOM) and NOM is more than an intercept:
    LL1 <- lapply(1:ncol(NOM), function(x) B1 * NOM[,x])
    B1 <- do.call(cbind, LL1)
    LL2 <- lapply(1:ncol(NOM), function(x) B2 * NOM[,x])
    B2 <- do.call(cbind, LL2)
  }
  ## update B1 and B2 with location effects (X):
  nbeta <- NCOL(X) - 1
  if(NCOL(X) > 1) {
    B1 <- cbind(B1, -X[, -1, drop = FALSE])
    B2 <- cbind(B2, -X[, -1, drop = FALSE])
  }
  dimnames(B1) <- NULL
  dimnames(B2) <- NULL
  list(B1=B1, B2=B2, o1=o1, o2=o2)
}

Deparse <-
  function(expr, width.cutoff = 500L, backtick = mode(expr) %in%
           c("call", "expression", "(", "function"),
           control = c("keepInteger", "showAttributes", "keepNA"),
           nlines = -1L)
### FIXME: test if formals(Deparse) == formals(deparse)??
  deparse(expr=expr, width.cutoff= width.cutoff, backtick=backtick,
          control=control, nlines=nlines)

