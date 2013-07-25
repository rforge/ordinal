clm <-
  function(formula, scale, nominal, data, weights, start, subset,
           doFit = TRUE, na.action, contrasts, model = TRUE,
           control = list(), ## tJac,
           link = c("logit", "probit", "cloglog", "loglog", "cauchit"),
           threshold = c("flexible", "symmetric", "symmetric2", "equidistant"), ...)
### deliberately no offset argument - include offset in the relevant
### formula/scale instead.
###
### FIXME: drop the tJac argument and let threshold accept it as a
### numeric matrix. Also test that ncol(tJac) <= nrow(tJac) or,
### perhaps better: rank(tJac) <= nlevels(y).
###
### FIXME: allow "threshold="fixed", theta=NULL" arguments to make it
### possible to fit certain kinds of models.
{
  ## Initial argument matching and testing:
  mc <- match.call(expand.dots = FALSE)
  link <- match.arg(link)
  threshold <- match.arg(threshold)
  ## check for presence of formula:
  if(missing(formula)) stop("Model needs a formula")
  if(missing(contrasts)) contrasts <- NULL
  ## set control parameters:
  control <- do.call(clm.control, c(control, list(...)))

  ## identify model as 'simple clm' or 'extended clm':
  Class <- if(!missing(scale) || link == "cauchit")
    c("eclm", "clm") else c("sclm", "clm")

  ## Compute: y, X, wts, off, mf:
  frames <- eclm.model.frame(mc, contrasts)

  ## Compute the transpose of the Jacobian for the threshold function,
  ## tJac and the names of the threshold parameters, alpha.names:
  ## if(missing(tJac)) tJac <- NULL
  ## if(!is.null(tJac) && missing(start))
  ##   stop("specify 'start' when supplying 'tJac'")
  ## frames$ths <- makeThresholds(frames$ylevels, threshold, tJac)
  frames$ths <- makeThresholds(frames$ylevels, threshold)

  ## Return model.frame?
  if(control$method == "model.frame") return(frames)

  ## Test column rank deficiency and possibly drop some
  ## parameters. Also set the lists aliased and coef.names:
  ## (X is with intercept at this point.)
  frames <- drop.cols(frames, drop.scale=FALSE, silent=TRUE)
### Note: intercept could be dropped from X and S in drop.cols?
### Currently they are dropped in eclm.newRho instead.

  ## Set envir rho with variables: B1, B2, o1, o2, wts, fitted:
  rho <- with(frames, {
      eclm.newRho(parent.frame(), y=y, X=X, NOM=frames$NOM, S=frames$S,
                  weights=wts, offset=off, S.offset=frames$S.off,
                  tJac=ths$tJac)
  })

  ## Set appropriate logLik and deriv functions in rho:
  if("eclm" %in% Class) {
    rho$clm.nll <- eclm.nll
    rho$clm.grad <- eclm.grad
    rho$clm.hess <- eclm.hess
  } else {
    rho$clm.nll <- clm.nll
    rho$clm.grad <- clm.grad
    rho$clm.hess <- clm.hess
  }

  ## Set starting values for the parameters:
  start <- set.start(rho, start=start, get.start=missing(start),
                     threshold=threshold, link=link, frames=frames)
  rho$par <- as.vector(start) ## remove attributes

  ## Set pfun, dfun and gfun in rho:
  setLinks(rho, link)

  ## Possibly return the environment rho without fitting:
  if(!doFit) return(rho)

  ## Fit the clm:
  if(control$method == "Newton") {
      fit <- clm.fit.NR(rho, control)
  } else
  fit <- clm.fit.optim(rho, control$method, control$ctrl)

  ## Modify and return results:
  res <- eclm.finalize(fit, weights=frames$wts,
                       coef.names=frames$coef.names,
                       aliased=frames$aliased)
  res$link <- link
  res$start <- start
  res$control <- control
  if(control$method == "Newton" &&
     !is.null(start.iter <- attr(start, "start.iter")))
    res$niter <- res$niter + start.iter
  res$threshold <- threshold
  res$call <- match.call()
  res$y.levels <- lev <- frames$ylevels
  res$tJac <- frames$ths$tJac
  rownames(res$tJac) <- paste(lev[-length(lev)], lev[-1], sep="|")
  colnames(res$tJac) <- frames$ths$alpha.names
  res$contrasts <- attr(frames$X, "contrasts")
  res$na.action <- attr(frames$mf, "na.action")
  res$terms <- frames$terms
  res$xlevels <- .getXlevels(res$terms, frames$mf)
  Theta.ok <- TRUE
  if(!is.null(frames$NOM)) {
      ## Save nominal information:
      res$nom.contrasts <- attr(frames$NOM, "contrasts")
      res$nom.terms <- frames$nom.terms
      res$nom.xlevels <- .getXlevels(res$nom.terms, frames$mf)
      ## get matrix of thresholds; Theta:
      Theta.list <-
          getThetamat(terms=res$nom.terms, alpha=res$alpha,
                      assign=attr(frames$NOM, "assign"),
                      contrasts=res$nom.contrasts, xlevels=res$nom.xlevels,
                      tJac=res$tJac)
### FIXME: cannot get Theta if some threshold parameters are aliased.
      ## Test that thresholds are increasing:
      if(all(is.finite(res$alpha))) {
          th.increasing <- apply(Theta.list$Theta, 1, function(th)
                                 all(diff(th) >= 0))
          if(!all(th.increasing))
              Theta.ok <- FALSE
      }
      res$Theta <- if(length(Theta.list) == 2)
          with(Theta.list, cbind(mf.basic, Theta)) else Theta.list$Theta
      res$alpha.mat <-
          matrix(res$alpha, ncol=ncol(res$tJac), byrow=TRUE)
      colnames(res$alpha.mat) <- colnames(res$tJac)
      rownames(res$alpha.mat) <- attr(frames$NOM, "orig.colnames")
  } else { ## no nominal effects:
      res$Theta <- res$alpha %*% t(res$tJac)
    ## colnames(res$Theta) <- names(res$alpha)
  }
  if(!is.null(frames$S)) {
    res$S.contrasts <- attr(frames$S, "contrasts")
    res$S.terms <- frames$S.terms
    res$S.xlevels <- .getXlevels(res$S.terms, frames$mf)
  }
  ## Check convergence:
### NOTE: we need to check convergence *after* generation of the Theta
### matrix, since here we check if thresholds are increasing.
  conv <- conv.check(res, Theta.ok=Theta.ok,
                     tol=control$tol)
  print.conv.check(conv, action=control$convergence) ## print convergence message
  res$vcov <- conv$vcov
  res$condHess <- conv$cond.H
  res$convergence <- conv[!names(conv) %in% c("vcov", "cond.H")]
  res$info <- with(res, {
      data.frame("link" = link,
                 "threshold" = threshold,
                 "nobs" = nobs,
                 "logLik" = formatC(logLik, digits=2, format="f"),
                 "AIC" = formatC(-2*logLik + 2*edf, digits=2,
                 format="f"),
                 "niter" = paste(niter[1], "(", niter[2], ")", sep=""),
### NOTE: iterations to get starting values for scale models *are*
### included here.
                 "max.grad" = formatC(maxGradient, digits=2,
                 format="e"),
                 "cond.H" = formatC(condHess, digits=1, format="e")
                 ## BIC is not part of output since it is not clear what
                 ## the no. observations are.
                 )
  })
  ## add model.frame to results list?
  if(model) res$model <- frames$mf

  res$par <- res$fitted <- res$niter <- NULL
  ## order elements of result alphabetically:
  res <- res[order(tolower(names(res)))]
  class(res) <- Class
  res
}

eclm.model.frame <- function(mc, contrasts) {
### mc - the matched call
### contrasts - contrasts for the model terms

    ## Collect all variables in a full formula:
  ## evaluate the formulae in the enviroment in which clm was called
  ## (parent.frame(2)) to get them evaluated properly:
  forms <- list(eval.parent(mc$formula, 2))
  if(!is.null(mc$scale)) forms$scale <- eval.parent(mc$scale, 2)
  if(!is.null(mc$nominal)) forms$nominal <- eval.parent(mc$nominal, 2)
  ## get the environment of the formula. If this does not have an
  ## enviroment (it could be a character), then use the parent frame.
  form.envir <-
    if(!is.null(env <- environment(forms[[1]]))) env
    else parent.frame(2)
  ## ensure formula, scale and nominal are formulas:
  ## forms <- lapply(forms, function(x) {
  ##   try(formula(deparse(x), env = form.envir), silent=TRUE) })
  for(i in 1:length(forms)) {
    forms[[i]] <- try(formula(deparse(forms[[i]]),
                              env = form.envir), silent=TRUE)
  }
  if(any(sapply(forms, function(f) class(f) == "try-error")))
    stop("unable to interpret 'formula', 'scale' or 'nominal'")
  ## collect all variables in a full formula:
  fullForm <- do.call("getFullForm", forms)
  ## set environment of 'fullForm' to the environment of 'formula':
  environment(fullForm) <- form.envir

  ## Extract the full model.frame(mf):
  m <- match(c("data", "subset", "weights", "na.action"),
             names(mc), 0)
  mf <- mc[c(1, m)]
  mf$formula <- fullForm
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  if(is.null(mf$data)) mf$data <- form.envir
  fullmf <- eval(mf, envir = parent.frame(2))
  mf$na.action <- "na.pass" ## filter NAs by hand below

  ## Extract X:
  ## get X from fullmf to handle e.g., NAs and subset correctly
  mf$formula <- forms[[1]]
  X.mf <- eval(mf, envir = parent.frame(2))
  X.terms <- attr(X.mf, "terms")
### FIXME: make sure that contrast for terms in X.terms are parsed to
### model.matrix here. Similar problem with S and NOM design matrices
### below.
  X <- model.matrix(X.terms, fullmf, contrasts)
  n <- nrow(X)
  ## Test for intercept in X:
  Xint <- match("(Intercept)", colnames(X), nomatch = 0)
  if(Xint <= 0) {
    X <- cbind("(Intercept)" = rep(1, n), X)
    warning("an intercept is needed and assumed in 'formula'",
            call.=FALSE)
  } ## intercept in X is guaranteed.
  off <- getOffset(X.mf)
  ## Filter NAs if any:
  if(!is.null(naa <- na.action(fullmf)))
    off <- off[-naa]

  wts <- getWeights(fullmf)
  ## Extract model response:
  y <- model.response(fullmf, "any") ## any storage mode
  if(!is.factor(y)) stop("response needs to be a factor", call.=FALSE)
  ## ylevels are the levels of y with positive weights
  ylevels <- levels(droplevels(y[wts > 0]))
  ## check that y has at least two levels:
  if(length(ylevels) == 1L)
      stop(gettextf("response has only 1 level ('%s'); expecting at least two levels",
                    ylevels), call.=FALSE)
  if(!length(ylevels))
      stop("response should be a factor with at least two levels")
  ## stop("response factor should have at least two levels")

  ## list of results:
  res <- list(y=y, ylevels=ylevels, X=X, wts=wts, off=off,
              mf=fullmf, terms=X.terms)

  ## Extract S (design matrix for the scale effects):
  if(!is.null(mc$scale)) {
    mf$formula <- forms$scale
    S.mf <- eval(mf, envir = parent.frame(2))
    if(!is.null(model.response(S.mf)))
      stop("response not allowed in 'scale'", call.=FALSE)
    res$S.terms <- attr(S.mf, "terms")
    S <- model.matrix(res$S.terms, fullmf, contrasts)
    ## Test for intercept in S:
    Sint <- match("(Intercept)", colnames(S), nomatch = 0)
    if(Sint <= 0) {
      S <- cbind("(Intercept)" = rep(1, n), S)
      warning("an intercept is needed and assumed in 'scale'",
              call.=FALSE)
    } ## intercept in S is guaranteed.
    res$S <- S
    res$S.off <- getOffset(S.mf)
    ## Filter NAs if any:
    if(!is.null(naa <- na.action(fullmf)))
      res$S.off <- res$S.off[-naa]
  }

  ## Extract NOM (design matrix for the nominal effects):
  if(!is.null(mc$nominal)) {
    mf$formula <- forms$nominal
    nom.mf <- eval(mf, envir = parent.frame(2))
    if(!is.null(model.response(nom.mf)))
      stop("response not allowed in 'nominal'", call.=FALSE)
    if(!is.null(model.offset(nom.mf)))
      stop("offset not allowed in 'nominal'", call.=FALSE)
    res$nom.terms <- attr(nom.mf, "terms")
    NOM <- model.matrix(res$nom.terms, fullmf, contrasts)
    NOMint <- match("(Intercept)", colnames(NOM), nomatch = 0)
    if(NOMint <= 0) {
      NOM <- cbind("(Intercept)" = rep(1, n), NOM)
      warning("an intercept is needed and assumed in 'nominal'",
              call.=FALSE)
    } ## intercept in NOM is guarantied.
    res$NOM <- NOM
  }

  ## return results:
  return(res)
  ## Note: X, S and NOM are with dimnames and intercepts are
  ## guaranteed. They may be column rank defecient.
}

eclm.newRho <-
  function(parent=parent.frame(), y, X, NOM=NULL, S=NULL, weights,
           offset, S.offset=NULL, tJac)
### Setting variables in rho: B1, B2, o1, o2, wts.
{
  rho <- new.env(parent = parent)
  ## Make B1, B2, o1, o2 based on y, X and tJac:
  ## rho <- list2env(getB(y=y, NOM=NOM, X=X, offset=offset,
  ## tJac=tJac), parent=parent)
  keep <- weights > 0
  y[!keep] <- NA
  y <- droplevels(y)
  ntheta <- nlevels(y) - 1
  y <- c(unclass(y))
  y[is.na(y)] <- 0
  n <- sum(keep)
  B2 <- 1 * (col(matrix(0, nrow(X), ntheta + 1)) == y)
  rho$o1 <- c(1e5 * B2[keep, ntheta + 1]) - offset[keep]
  rho$o2 <- c(-1e5 * B2[keep, 1]) - offset[keep]
  B1 <- B2[keep, -(ntheta + 1), drop = FALSE]
  B2 <- B2[keep, -1, drop = FALSE]
  ## adjust B1 and B2 for structured thresholds:
  rho$B1 <- B1 %*% tJac
  rho$B2 <- B2 %*% tJac
  ## update B1 and B2 with nominal effects:
  if(NCOL(NOM) > 1) { ## !is.null(NOM) && ncol(NOM) > 1) {
    ## if !is.null(NOM) and NOM is more than an intercept:
    LL1 <- lapply(1:ncol(NOM), function(x) rho$B1 * NOM[keep, x])
    rho$B1 <- do.call(cbind, LL1)
    LL2 <- lapply(1:ncol(NOM), function(x) rho$B2 * NOM[keep, x])
    rho$B2 <- do.call(cbind, LL2)
  }
  ## update B1 and B2 with location effects (X):
  nbeta <- NCOL(X) - 1
  if(nbeta > 0) {
    rho$B1 <- cbind(rho$B1, -X[keep, -1, drop = FALSE])
    rho$B2 <- cbind(rho$B2, -X[keep, -1, drop = FALSE])
  }
  dimnames(rho$B1) <- NULL
  dimnames(rho$B2) <- NULL
  rho$n.psi <- ncol(rho$B1) ## no. linear model parameters
  rho$k <- 0
  ## there may be scale offset without scale predictors:
  rho$sigma <- rho$Soff <-
    if(is.null(S.offset)) rep(1, n) else exp(S.offset[keep])
  ## save scale model:
  if(!is.null(S)) {
    rho$S <- S[keep, -1, drop=FALSE]
    dimnames(rho$S) <- NULL
    rho$k <- ncol(rho$S) ## no. scale parameters
  }
  rho$has.scale <- ## TRUE if scale has to be considered.
    (!is.null(S) || any(S.offset != 0))
  ## initialize fitted values and weights:
  rho$fitted <- numeric(length = n)
  rho$wts <- weights[keep]
  ## return:
  return(rho)
}

clm.fit.optim <-
  function(rho, method = c("ucminf", "nlminb", "optim"), control=list())
{
  method <- match.arg(method)
  ## optimize the likelihood:
  optRes <-
    switch(method,
           "nlminb" = nlminb(rho$par,
             function(par) eclm.nll(rho, par),
             function(par) eclm.grad2(rho, par),
             control=control),
           "ucminf" = ucminf(rho$par,
             function(par) eclm.nll(rho, par),
             function(par) eclm.grad2(rho, par),
             control=control),
           "optim" = optim(rho$par,
             function(par) eclm.nll(rho, par),
             function(par) eclm.grad2(rho, par),
             method="BFGS",
             control=control),
           )
  ## save results:
  rho$par <- optRes[[1]]
  res <- list(par = rho$par,
              logLik = -eclm.nll(rho),
              gradient = eclm.grad(rho),
              Hessian = eclm.hess(rho),
              fitted = rho$fitted)
  res$maxGradient = max(abs(res$gradient))
  res$optRes <- optRes
  res$niter <- switch(method, "nlminb" = optRes$evaluations,
                      "ucminf" = c(optRes$info["neval"], 0),
                      "optim" = optRes$counts)
  res$convergence <-
    switch(method, "nlminb" = optRes$convergence,
           "ucminf" = optRes$convergence,
           "optim" = optRes$convergence)

  return(res)
}

eclm.nll <- function(rho, par) {
  if(!missing(par)) rho$par <- par
  with(rho, {
      if(k > 0)
      sigma <- Soff * exp(drop(S %*% par[n.psi + 1:k]))
### NOTE: we have to divide by sigma even if k=0 since there may be an
### offset but no predictors in the scale model:
    eta1 <- (drop(B1 %*% par[1:n.psi]) + o1)/sigma
    eta2 <- (drop(B2 %*% par[1:n.psi]) + o2)/sigma
  })
### NOTE: getFitted is not found from within rho, so we have to
### evalueate it outside of rho
  rho$fitted <- getFittedC(rho$eta1, rho$eta2, rho$link)
  if(all(is.finite(rho$fitted)) && all(rho$fitted > 0))
### NOTE: Need test here because some fitted <= 0 if thresholds are
### not ordered increasingly.
    -sum(rho$wts * log(rho$fitted))
  else Inf
}

eclm.grad <- function(rho) {
### requires that eclm.nll has been called prior to
### eclm.grad.
  with(rho, {
    p1 <- dfun(eta1)
    p2 <- dfun(eta2)
    wtpr <- wts/fitted
    C2 <- B1*p1/sigma - B2*p2/sigma
    if(k <= 0) return(-crossprod(C2, wtpr))
    C3 <- -(eta1 * p1 - eta2 * p2) * S
    return(-crossprod(cbind(C2, C3), wtpr))
### NOTE: C2 and C3 are used by eclm.hess
  })
}

eclm.grad2 <- function(rho, par) {
### does not require that eclm.nll has been called prior to
### eclm.grad.
  eclm.nll(rho, par)
  eclm.grad(rho)
}

eclm.hess <- function(rho) {
### requires that eclm.grad has been called prior to this.
  with(rho, {
    g1 <- gfun(eta1)
    g2 <- gfun(eta2)
    wtprpr <- wtpr/fitted ## Phi3
    dg.psi <- crossprod(B1 * gfun(eta1) * wtpr / sigma^2, B1) -
      crossprod(B2 * gfun(eta2) * wtpr / sigma^2, B2)
    ## upper left:
    D <- dg.psi - crossprod(C2, (C2 * wtprpr))
    if(k <= 0) return(-D) ## no scale predictors
    ## upper right (lower left transpose):
    wtprsig <- wtpr/sigma
    epg1 <- p1 + g1*eta1
    epg2 <- p2 + g2*eta2
    Et <- crossprod(B1, -wtprsig * epg1 * S) -
      crossprod(B2, -wtprsig * epg2 * S) -
        crossprod(C2, wtprpr * C3)
    ## lower right:
    F <- -crossprod(S, wtpr * ((eta1*p1 - eta2*p2)^2 / fitted -
                               (eta1*epg1 - eta2*epg2)) * S)
    ## combine and return hessian:
    H <- rbind(cbind(D    , Et),
               cbind(t(Et), F))
    return(-H)
  })
}


##################################################################
### NEW VERSION OF FITTING ALGORITHM:
##################################################################

clm.fit.NR <-
  function(rho, control = list())
### The main work horse: Where the actual fitting of the clm goes on.
### Fitting the clm via modified Newton-Raphson with step halving.

### -------- Assumes the existence of the following functions:
### eclm.nll - negative log-likelihood
### eclm.grad - gradient of nll wrt. par
### eclm.hess - hessian of nll wrt. par
### Trace - for trace information
{
    control <- do.call(clm.control, control)
    stepFactor <- 1
    innerIter <- modif.iter <- 0L
    conv <- 2L  ## Convergence flag (iteration limit reached)
    nll <- rho$clm.nll(rho)
    if(!is.finite(nll))
        stop("Non-finite log-likelihood at starting value")
    ## do.newton <-
    ## rel.conv <- FALSE
    ## stephalf <- TRUE

    ## Newton-Raphson algorithm:
    for(i in 1:(control$maxIter + 1L)) {
        gradient <- rho$clm.grad(rho)
        maxGrad <- max(abs(gradient))
        if(control$trace > 0) {
            Trace(iter=i+innerIter-1, stepFactor, nll, maxGrad,
                  rho$par, first=(i==1))
            if(control$trace > 1 && i > 1) {
                cat("\tgrad: ")
                cat(paste(formatC(gradient, digits=3, format="e")))
                cat("\n\tstep: ")
                cat(paste(formatC(-step, digits=3, format="e")))
                cat("\n\teigen: ")
                cat(paste(formatC(eigen(hessian, symmetric=TRUE,
                                        only.values=TRUE)$values, digits=3,
                                  format="e")))
                cat("\n")
            }
        }
        abs.conv <- (maxGrad < control$gradTol)
        hessian <- rho$clm.hess(rho)
        ## Compute cholesky factor of Hessian: ch = Ut U
        ch <- try(chol(hessian), silent=TRUE)
### NOTE: solve(hessian, gradient) is not good enough because it will
### compute step for negative-definite Hessians and we don't want
### that.

### FIXME: What if Hessian is closely singular but slightly positive?
### Could we do something better in that case?
        if(inherits(ch, "try-error")) {
            if(abs.conv) { ## step.ok not true.
                conv <- 1L
                break ## cannot meet relative criterion.
            }
            ## If Hessian is non-positive definite:
            min.ev <- min(eigen(hessian, symmetric=TRUE,
                                only.values=TRUE)$values)
            inflation.factor <- 1
            ## Inflate diagonal of Hessian to make it positive definite:
            inflate <- abs(min.ev) + inflation.factor
            hessian <- hessian + diag(inflate, nrow(hessian))
            if(control$trace > 0)
                cat(paste("Hessian is singular at iteration", i-1, "inflating diagonal with",
                          formatC(inflate, digits=5, format="f"), "\n"))
            ch <- try(chol(hessian), silent=TRUE)
            if(inherits(ch, "try-error"))
                stop(gettextf("Cannot compute Newton step at iteration %d",
                              i-1), call.=FALSE)
            modif.iter <- modif.iter + 1L
            ## do.newton <- FALSE
        } else
            modif.iter <- 0L
        if(modif.iter >= control$maxModIter) {
            conv <- 4L
            break
        }

        ## solve U'y = g for y, then
        ## solve U step = y for step:
        step <- c(backsolve(ch, backsolve(ch, gradient, transpose=TRUE)))
        rel.conv <- (max(abs(step)) < control$relTol)
        ## Test if step is in a descent direction -
        ## otherwise use step <- grad / max|grad|:
        ## if(crossprod(gradient, step) < 0) {
        ##     if(control$trace > 0)
        ##         cat("Newton step is not in descent direction; using gradient instead\n")
        ##     step <- c(gradient / max(abs(gradient)))
        ## } else

        ## update parameters:
        rho$par <- rho$par - stepFactor * step
        nllTry <- rho$clm.nll(rho)
        lineIter <- 0
        stephalf <- (nllTry > nll)
### NOTE: sometimes nllTry > nll just due to noise, so we also check
### reduction in gradient for small diffs:
        if(stephalf && abs(nll - nllTry) < 1e-10)
            stephalf <- maxGrad < max(abs(rho$clm.grad(rho)))

        ## Assess convergence:
        ## (only attempt to sattisfy rel.conv if abs.conv is true and
        ## it is possible to take the full newton step)
        ## if(abs.conv && !step.ok) {
        if(abs.conv && stephalf) {
            conv <- 1L
            ## we need to step back to the par for which abs.conv
            ## was true:
            rho$par <- rho$par + stepFactor * step
            rho$clm.nll(rho)
            break
        }
        if(abs.conv && rel.conv) {
            conv <- 0L
            ## no need to step back as stephalf was false so the new
            ## par are just better.
            break
        }
        ## Step halving if nll increases:
        while(stephalf) {
            stepFactor <- stepFactor/2
            rho$par <- rho$par + stepFactor * step
            nllTry <- rho$clm.nll(rho)
            lineIter <- lineIter + 1
            if(control$trace > 0) {
                cat("step halving:\n")
                cat("nll reduction: ", formatC(nll - nllTry, digits=5, format="e"), "\n")
                Trace(i+innerIter-1, stepFactor, nll, maxGrad,
                      rho$par, first = FALSE)
            }
            if(lineIter > control$maxLineIter){
                conv <- 3L
                break
            }
            innerIter <- innerIter + 1
            stephalf <- (nllTry > nll)
            if(stephalf && abs(nll - nllTry) < 1e-10)
                stephalf <- (maxGrad < max(abs(rho$clm.grad(rho))))
        } ## end step halving
        if(conv == 3L) break

        if(control$trace > 0)
            cat("nll reduction: ", formatC(nll - nllTry, digits=5, format="e"), "\n")
        nll <- nllTry
        ## Double stepFactor if needed:
        stepFactor <- min(1, 2 * stepFactor)
    } ## end Newton iterations
    message <- switch(as.character(conv),
                      "0" = "Absolute and relative convergence criteria were met",
                      "1" = "Absolute convergence criterion was met, but relative criterion was not met",
                      "2" = "iteration limit reached",
                      "3" = "step factor reduced below minimum",
                      "4" = "maximum number of consecutive Newton modifications reached")
    if(conv <= 1L && control$trace > 0) {
        cat("\nOptimizer converged! ", message, fill = TRUE)
    }
    if(conv > 1 && control$trace > 0) {
        cat("\nOptimization failed ", message, fill = TRUE)
    }
    ## return results:
    res <- list(par = rho$par,
                gradient = c(rho$clm.grad(rho)), ##as.vector(gradient),
                ## Hessian = hessian,
                Hessian = rho$clm.hess(rho), ## ensure hessian is evaluated
                ## at optimum
                logLik = -nll,
                convergence = conv,
                ## 0: abs and rel criteria meet
                ## 1: abs criteria meet, rel criteria not meet
                ## 2: iteration limit reached
                ## 3: step factor reduced below minium
                message = message,
                maxGradient = maxGrad,
                niter = c(outer = i-1, inner = innerIter),
                fitted = rho$fitted)
    return(res)
}

