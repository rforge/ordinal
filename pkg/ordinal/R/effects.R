gof2 <-
  function(object, test = c("Pearson", "Deviance"), ...)
{
  ## FIXME: Since PRODID is contained in PROD, should have n.PRODID
  ## rows and not n.PRODID * n.PROD rows.
  xlevels <- object$xlevels
  contrasts <- object$contrasts
  Terms <- delete.response(terms(object))
  ## FIXME: include terms/variables from scale and nominal parts here as
  ## well. 
  factor.table <- attr(Terms, "factors")
  all.varnm <- rownames(factor.table)
  var.classes <- attr(Terms, "dataClasses")[all.varnm]
  numeric.var <- which(var.classes != "factor")
  ## check for numeric.var in all three formulae:
  if(length(numeric.var))
    stop("cannot cope with numeric variables")
  factor.var <- which(var.classes == "factor")
  if(!length(factor.var))
    stop("need factor variables in model")
  
  factor.varnm <- names(var.classes)[factor.var]
  ## need union of all xlevels here:
  xlev <- xlevels[factor.varnm]
  ## minimal complete model frame:
  (mf.basic <- do.call(expand.grid, xlev))
  ## minimal complete design matrix:
  prob <- predict.clm(object, newdata=mf.basic)$fit
  prob.table <- cbind(mf.basic, prob)
  
  ## Getting observed frequencies:
  mf <- model.frame(object)
  wts <- if("(weights)" %in% names(mf))
    mf[["(weights)"]] else  rep(1, nrow(mf)) 
  sp.wts <- split(wts, mf[factor.varnm])
  sp.y <- split(mf[[response.name(terms(object))]], mf[factor.varnm])
  obs.tabs <- lapply(1:length(sp.y), function(i) {
    tab <- tapply(sp.wts[[i]], sp.y[[i]], sum)
    tab[is.na(tab)] <- 0
    tab
  })
  observed <- do.call(rbind, obs.tabs)
  data.table <- cbind(mf.basic, as.data.frame(observed))
  ## Getting expected frequencies:
  margin <- rowSums(observed)
  expected <-
    lapply(1:length(margin), function(i) prob[i, ] * margin[i])
  expected <- do.call(rbind, expected)
  if(any(expected < 5))
    warning("Chisq approximation may be inaccurate")
  ## computing gof test:
  test <- match.arg(test)
  if(test == "Pearson") {
    resid <- (observed - expected) / sqrt(expected)
    X2 <- sum(resid^2)
  }
  if(test == "Deviance"){
    resid <- NULL
    nz <- observed > 0
    X2 <- 2 * sum(observed[nz] * log(observed[nz]/expected[nz]))
  }
  resdf <- nrow(observed) * (ncol(observed) - 1) - object$edf
  P <- pchisq(X2, df=resdf, lower.tail=FALSE)
  ## results:
  res <- list(data.table=mf.basic, prob.table=prob,
              observed=observed, expected=expected, test=test,
              statistic=X2, p.value=P, residual.df=resdf, resid=resid)  
  res  
}

getTheta <-
  function(object, type=c("Theta", "Alpha"), include.cell.id=TRUE,
           ...) 
### Predict thresholds for all (combination of levels) of categorical
### nominal variables. Continuous variables are ignored. 
{
  type <- match.arg(type)
  if(is.null(object$call$nominal))
    return(object$alpha)
  Terms <- object$nom.terms
  xlev <- object$nom.xlevels
  con <- object$nom.contrasts
  mf.basic <- do.call(expand.grid, xlev)
  X <- model.matrix(Terms, data=mf.basic, contrasts=con)
  Theta <- matrix(object$alpha, ncol=ncol(object$tJac), byrow=TRUE)
  Theta <- apply(Theta, 2, function(th) X %*% th)
  if(type == "Theta") {
    Theta <- t(apply(Theta, 1, function(th) object$tJac %*% th))
    lev <- object$y.levels
    colnames(Theta) <- rownames(object$tJac)
  } else {
    colnames(Theta) <- colnames(object$tJac)
  }
  Theta <- as.data.frame(Theta)
  if(include.cell.id)
    Theta <- cbind(mf.basic, Theta)
  Theta
}

getalpha.names <- function(clm_object, ...)
{
  obj <- clm_object
  obj$tJac
}
