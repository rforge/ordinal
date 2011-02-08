### Make small self-contained functions.
### Keep code in only one place - no code should be dublicated.
### This probably means that all (most) functions should be globally
### visible and not hidden inside functions.

### Rule: a function taking an environment as argument will update
### some of the elements of the environment, but not introduce new
### ones (without removing them again). New values are returned from
### the function instead.

### FILES:
### Main fcts:
### links
### start
### threshold(?)

##  start.type = 1L: only thresholds
##  2L: thresholds and beta = 0
##  3l: thresholds and beta via glm.fit
### Need y in rho to set starting values via binomial glm.


### -------- Notes ------------
### - change 'scale' to 'log scale' in print and summary
### - drop observations with zero or negative weights - at least give
###   a warning.
### - alter print.summary.clmm to display data information as in
###   print.lmer
### - implement fixef() and ranef()
### - use ginv in NR alg to accomodate non-identifiability?
### - be more particular about method dispatch and class inheritance
### - to avoid e.g. clm.dropterm being called by a clmm object.
### - test rank of B1 / B2 when some threshold parameters are
###   undefined? Is this enough of to ensure identifiability (by
###   dropping some coef)?
### - format cond(Hess) in scientific format - it is the order of
###   magnitude that is relevant.

clm <-
  function(formula, data, weights, start, subset, doFit = TRUE,
           na.action, contrasts, model = TRUE, control,
           link = c("logit", "probit", "cloglog", "loglog", "cauchit"), 
           threshold = c("flexible", "symmetric", "equidistant"), ...)
{
  ## Initial argument matching and testing:
  mc <- match.call(expand.dots = FALSE)
  link <- match.arg(link)
  threshold <- match.arg(threshold)
  ## check for presence of formula:
  if(missing(formula)) stop("Model needs a formula")
  if(missing(contrasts)) contrasts <- NULL
  
  ## set control parameters:
  if(missing(control)) control <- clm.control(...)
  if(!setequal(names(control), c("method", "trace", "maxIter",
                                 "gradTol", "maxLineIter")))
    stop("specify 'control' via clm.control()")

  ## make data a data.frame:
  ##    if(missing(data)) mc$data <- environment(formula)
  ##    if(is.matrix(eval.parent(mc$data)))
  ##      mc$data <- as.data.frame(mc$data)
  
  ## Compute: y, X, wts, off, mf:
  frames <- clm.model.frame(mc, contrasts)
  if(control$method == "model.frame") return(frames)

  ## Test rank deficiency and possibly drop some parameters:
  ## X is guarantied to have an intercept at this point.
  frames$X <- dropCoef(frames$X)
  
  ## Compute the transpose of the Jacobian for the threshold function,
  ## tJac and the names of the threshold parameters, alphaNames:
  ths <- makeThresholds(frames$y, threshold)

  ## set envir; rho with variables: B1, B2, o1, o2, wts, pfun,
  ## dfun, gfun, fitted, control/ctrl:
  rho <- with(frames, {
    clm.newRho(parent.frame(), y = y, X = X, weights = wts,
               offset = off, link = link,  tJac = ths$tJac) })

  ## set starting values for the parameters:
  if(missing(start))
    start <- clm.start(frames$y, frames$X, has.intercept = TRUE,
                       threshold = threshold)
  stopifnot(is.numeric(start) && 
            length(start) == (ths$nalpha + ncol(frames$X) - 1) )
  rho$par <- start
### FIXME: if(link == "cauchit") fit probit clm and use coef as
### start. 

  ## possibly return the environment, rho without fitting:
  if(!doFit) return(rho)

  ## fit the clm:
  control$method <- NULL
  fit <- clm.fit.env(rho, control)
### FIXME: add arg non.conv = c("error", "warn", "message") ?

  ## Modify and return results:
  res <- clm.finalize(fit, weights = frames$wts, mc = match.call(),
                      contrasts = contrasts,
                      tJac = ths$tJac,
                      na.action = attr(frames$mf, "na.action"),
                      terms = frames$terms, 
                      alphaNames = ths$alphaNames,
                      betaNames = colnames(frames$X)[-1])

  ## add model.frame to results list?
  if(model) res$model <- frames$mf
  
  return(res)
}

clm.finalize <- function(fit, weights, mc, contrasts, tJac, na.action, 
                         terms, alphaNames, betaNames)
{   
  nalpha <- length(alphaNames)
  nbeta <- length(betaNames)
  stopifnot(length(fit$par) == nalpha + nbeta)
  
  fit <- within(fit, {
    alpha <- par[1:nalpha]
    names(alpha) <- alphaNames
    beta <- if(nbeta > 0) par[nalpha + 1:nbeta]
    else numeric(0)
    names(beta) <- betaNames
    coefficients <- c(alpha, beta)
    
    names(gradient) <- names(coefficients)
    dimnames(Hessian) <- list(names(coefficients),
                              names(coefficients))
    edf <- length(coefficients) ## estimated degrees of freedom
    nobs <- sum(weights)
    n <- length(weights)
    fitted.values <- fitted
    df.residual = nobs - edf

    tJac <- tJac
    terms <- terms
    contrasts <- contrasts
    na.action <- na.action
    call <- mc
    logLik <- -fit$nll

    rm(list = c("par"))
  })
  class(fit) <- "clm"
  return(fit)
}

dropCoef <- function(X) {
### For a clm it is assumed that X has an intercept
  stopifnot(is.matrix(X))
  
  qr.X <- qr(X, tol = 1e-7, LAPACK = FALSE)
  if(qr.X$rank == ncol(X)) return(X)

  warning(gettextf("design is rank deficient so dropping %d coef",
                   ncol(X) - qr.X$rank), call. = FALSE)
  ## warning(gettextf("design is rank deficient so dropping some coef"))
  ## take the columns correponding to the first qr.x$rank pivot
  ## elements of X:
  newX <- X[, qr.X$pivot[1:qr.X$rank], drop = FALSE]
  if(qr.X$rank != qr(newX)$rank)
    stop(gettextf("determination of full rank design matrix failed"),
         call. = FALSE)
### Should the rank deficiency be testet on B1 and B2?
  return(newX) 
}

clm.model.frame <- function(mc, contrasts) {
### mc - the matched call
### contrasts - contrasts for the fixed model terms

  ## Extract model.frame(mf):
  m <- match(c("formula", "data", "subset", "weights", "na.action",
               "offset"), names(mc), 0)
  mf <- mc[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, envir = parent.frame(2))

  ## Extract model response:
  y <- model.response(mf, "any")
### FIXME: why "any"?
  if(!is.factor(y))
    stop("response needs to be a factor")

  ## Extract X:
  terms <- attr(mf, "terms")
  X <- model.matrix(terms, mf, contrasts)
  n <- nrow(X)
  ## Test for intercept in X:
  Xint <- match("(Intercept)", colnames(X), nomatch = 0)
  if(Xint <= 0) {
    X <- cbind("(Intercept)" = rep(1, n), X)
    warning("an intercept is needed and assumed")
  } ## intercept in X is guarantied.

  ## Extract the weights and offset:
  if(is.null(wts <- model.weights(mf))) wts <- rep(1, n)
  if(is.null(off <- model.offset(mf))) off <- rep(0, n)
  
  ## check weights and offset:
  if (any(wts <= 0))
    stop(gettextf("non-positive weights are not allowed"))
### FIXME: Wouldn't it be usefull to just remove any observations with
### zero weights?
  if(length(wts) && length(wts) != NROW(y))
    stop(gettextf("number of weights is %d should equal %d
(number of observations)", length(wts), NROW(y)))
  if(length(off) && length(off) != NROW(y))
    stop(gettextf("number of offsets is %d should equal %d
(number of observations)", length(off), NROW(y)))
  
  ## return list of results:
  res <- list(y = y, X = X, wts = as.double(wts),
              off = as.double(off), mf = mf, terms = terms)
  ## Note: X is with dimnames and an intercept is guarantied.
  return(res)
}

setLinks <- function(rho, link) {
### The Aranda-Ordaz and log-gamma links are not supported in this
### version of clm.
  rho$pfun <- switch(link,
                     logit = plogis,
                     probit = pnorm,
                     cloglog = pgumbel,
                     cauchit = pcauchy,
                     loglog = pgumbel2,
                     "Aranda-Ordaz" = function(x, lambda) pAO(x, lambda),
                     "log-gamma" = function(x, lambda) plgamma(x, lambda))
  rho$dfun <- switch(link,
                     logit = dlogis,
                     probit = dnorm,
                     cloglog = dgumbel,
                     cauchit = dcauchy,
                     loglog = dgumbel2,
                     "Aranda-Ordaz" = function(x, lambda) dAO(x, lambda),
                     "log-gamma" = function(x, lambda) dlgamma(x, lambda))
  rho$gfun <- switch(link,
                     logit = glogis,
                     probit = function(x) -x * dnorm(x),
                     loglog = function(x) -ggumbel(-x),
                     cloglog = ggumbel,
                     cauchit = gcauchy,
                     "Aranda-Ordaz" = function(x, lambda) gAO(x, lambda), ## shouldn't happen
                     "log-gamma" = function(x, lambda) glgamma(x, lambda)
                     )
  rho$link <- link
}

start.threshold <-
  function(y, threshold = c("flexible", "symmetric", "equidistant")) 
### args:
### y - model response, a factor with at least two levels
### threshold - threshold structure, character.
{
  ## match and test arguments:
  threshold <- match.arg(threshold)
  stopifnot(is.factor(y) && nlevels(y) >= 2)
  ntheta <- nlevels(y) - 1L
  if(threshold %in% c("symmetric", "equidistant") && nlevels(y) < 3)
    stop(gettextf("symmetric and equidistant thresholds are only
meaningful for responses with 3 or more levels"))
  
  ## default starting values:
  start <- qlogis((1:ntheta) / (ntheta + 1) ) # just a guess
  
  ## adjusting for threshold functions:
  if(threshold == "symmetric" && ntheta %% 2) { ## ntheta odd >= 3
    nalpha <- (ntheta + 1) / 2
    start <- c(start[nalpha], diff(start[nalpha:ntheta])) ## works for
    ## ntheta >= 1
  }
  if(threshold == "symmetric" && !ntheta %% 2) {## ntheta even >= 4
    nalpha <- (ntheta + 2) / 2
    start <- c(start[c(nalpha - 1, nalpha)],
               diff(start[nalpha:ntheta])) ## works for ntheta >= 2
  }
  if(threshold == "equidistant")
    start <- c(start[1], mean(diff(start)))

  ## return starting values for the threshold parameters:
  return(as.vector(start))
}

start.beta <- function(X, has.intercept = TRUE)
  return(rep(0, NCOL(X) - has.intercept))

clm.start <- function(y, threshold, X, has.intercept = TRUE)
  return(c(start.threshold(y, threshold),
           start.beta(X, has.intercept)))  

### Will these very simple starting values be enough for the NR alg to
### converge robustly to the optimum?  


### Parameters in rho needed to optimize a clm:
### Initial:
### B1, B2, nxi, p, par(initial), o1, o2, wts, pfun, dfun, gfun,
### control 

### Generated during fitting:
### eta1, eta2, pr, wtpr, p1, p2, dS.psi, dpi.psi

### Variables needed to set starting values for a clm:
### Thresholds: ntheta/nlev.y/y, threshold
### Regression par: y, X, wts, off, link, OR regParStart <- rep(0, p)

clm.newRho <-
  function(parent, y, X, weights, offset, link, tJac)
### Setting variables in rho: B1, B2, o1, o2, wts, pfun,
### dfun, gfun
{
  rho <- new.env(parent = parent)
  
  ## Set pfun, dfun and gfun in rho:
  setLinks(rho, link)
### FIXME: Should links be set like tJac etc.

  ## Make B1, B2, o1, o2 based on y, X and tJac:
  ntheta <- nlevels(y) - 1
  n <- nrow(X)
  B2 <- 1 * (col(matrix(0, n, ntheta + 1)) == c(unclass(y)))
  rho$o1 <- c(100 * B2[, ntheta + 1]) - offset
  rho$o2 <- c(-100 * B2[,1]) - offset
  B1 <- B2[,-(ntheta + 1), drop=FALSE]
  B2 <- B2[,-1, drop=FALSE]
  rho$B1 <- B1 %*% tJac
  rho$B2 <- B2 %*% tJac
  nbeta <- NCOL(X) - 1
  if(nbeta > 0) {
    ## Update B1 and B2 with X minus intercept:
    rho$B1 <- cbind(rho$B1, -X[, -1, drop = FALSE])
    rho$B2 <- cbind(rho$B2, -X[, -1, drop = FALSE])
  }
  dimnames(rho$B1) <- NULL
  dimnames(rho$B2) <- NULL
  
  rho$fitted <- numeric(length = n)
  rho$wts <- weights

  return(rho)
}

clm.fit <-
  function(y, X, weights = rep(1, nrow(X)), offset = rep(0, nrow(X)),
           control = list(), start, 
           link = c("logit", "probit", "cloglog", "loglog"), 
           threshold = c("flexible", "symmetric", "equidistant"))
{
  ## Initial argument matching and testing:
  threshold <- match.arg(threshold)
  link <- match.arg(link)
  control <- do.call(clm.control, control)
  stopifnot(is.factor(y) &&
            is.matrix(X))
  stopifnot(length(y) == nrow(X) &&
            length(y) == length(weights) &&
            length(y) == length(offset))
  ## Test for intercept in X:
  Xint <- match("(Intercept)", colnames(X), nomatch = 0)
  if(Xint <= 0) {
    X <- cbind("(Intercept)" = rep(1, nrow(X)), X)
    warning("an intercept is needed and assumed")
  } ## intercept in X is guarantied.

  ## ensure X has full rank, generate threshold structure and set the
  ## rho environment:
  X <- dropCoef(X)
  ths <- makeThresholds(y, threshold)
  rho <- clm.newRho(parent.frame(), y = y, X = X, weights = weights,
                    offset = offset, link = link, tJac = ths$tJac)

  ## set starting values for the clm:
  if(missing(start))
    start <- clm.start(y, X, has.intercept = TRUE, threshold = threshold)
  stopifnot(is.numeric(start) && 
            length(start) == (ths$nalpha + ncol(X) - 1) )
  rho$par <- start
  
  ## fit the clm:
  fit <- clm.fit.env(rho, control = control)
  return(fit)
### FIXME: should 'par' be 'coefficients'?
### FIXME: should fit contain 'vcov'?
}
  
clm.fit.env <-
  function(rho, control = list(trace = 0, maxIter = 100,
                  gradTol = 1e-4, maxLineIter = 10))
### The main work horse: Where the actual fitting of the clm goes on.
### Fitting the clm via Newton-Raphson with step halfing. 

### -------- Assumes the existence of the following functions:
### clm.nll - negative log-likelihood
### clm.grad - gradient of nll wrt. par
### clm.hess - hessian of nll wrt. par
### Trace - for trace information

### FIXME: use do.call(clm.control, control) to ensure none are
### missing? 
{
  stepFactor <- 1
  innerIter <- 0
  conv <- 1  ## Convergence flag
  message <- "iteration limit reached"
  nll <- clm.nll(rho)
  if(!is.finite(nll))
    stop("Non-finite log-likelihood at starting value")
  gradient <- clm.grad(rho)
  maxGrad <- max(abs(gradient))
  hessian <- clm.hess(rho)
  if(control$trace > 0)
    Trace(iter=0, stepFactor, nll, maxGrad, rho$par, first=TRUE)
  
  ## Newton-Raphson algorithm:
  for(i in 1:control$maxIter) {
    if(maxGrad < control$gradTol) {
      message <- "max|gradient| < tol, so current iterate is probably solution"
      if(control$trace > 0)
        cat("\nOptimizer converged! ", "max|grad|:",
            maxGrad, message, fill = TRUE)
      conv <- 0
      break
    } ## end convergence test

    ## Compute Newton step and update parameters
    step <- .Call("La_dgesv", hessian, gradient, .Machine$double.eps,
                  PACKAGE = "base") ## solve H*step = g for 'step'
### FIXME: possibly use ginv instead for unidentifiable model?
    step <- as.vector(step)
    rho$par <- rho$par - stepFactor * step
    nllTry <- clm.nll(rho)
    lineIter <- 0

    ## Step halfing if nll increases:
    while(nllTry > nll) {
      ## Step-half if nll increases or if nll is flat and maxGrad is
      ## not decreased. Thus if nll is flat, but maxGrad is decreased,
      ## the full step is justified:
      ## while(nllTry > nll || nllTry == nll && maxGradTry >= maxGrad) {
### Can nll increase while maxGrad decreases?

### It seems peculiar to ponder about a change in parameters that is
### so small that the likelihood is entirely flat. Still, the
### parameters can be determined with higher precision than what is
### visible from changes in nll.

### FIXME: Takes the step even if nllTry == nll - is this reasonable?
### This does happen with the rhyme data set. Apparently this is due
### to discreteness in the log-likelhood function, so even though a
### step can be taken that reduces max|grad|, the nll is flat and not
### reduced by the step.
### Could and should this be avoided by terminating if the step is not
### large enough? e.g., max(abs(step)) > 1e-8
### In a very well determined problem, the step can be very small -
### not if the gradient is very small, but if the Hessian is very
### large. Effectively par are at the MLEs.
      ## if(nllTry == nll) {
      ##   gradientTry <- clm.grad(rho)
      ##   if(maxGrad <= max(abs(clm.grad(rho)))) {
      ##     conv <- 3
      ##     message <- ""
      ## }
      stepFactor <- stepFactor/2
      rho$par <- rho$par + stepFactor * step
      nllTry <- clm.nll(rho)
      lineIter <- lineIter + 1
      if(control$trace > 0)
        Trace(i+innerIter, stepFactor, nll, maxGrad,
              rho$par, first = FALSE)
      if(lineIter > control$maxLineIter){
        message <- "step factor reduced below minimum"
        conv <- 2
        break
      }
      innerIter <- innerIter + 1
    } ## end step halfing
    if(conv == 2) break

    ## Update nll, gradient, maxGrad and hessian:
    gradient <- clm.grad(rho)
    maxGrad <- max(abs(gradient))
    hessian <- clm.hess(rho)
    nll <- nllTry
    if(control$trace > 0)
      Trace(iter=i+innerIter, stepFactor, nll,
            maxGrad, rho$par, first = FALSE)
    ## Double stepFactor if needed:
    stepFactor <- min(1, 2 * stepFactor)
  } ## end Newton iterations

### FIXME: warn if fitted values are close to 0 or 1
  
  if(conv > 0) { ## optimization failed
    if(control$trace > 0) cat(message, fill = TRUE)
    stop(gettextf("optimization failed: %s", message), call. = FALSE)
  }
  
  ## return results:
  res <- list(par = rho$par,
              gradient = as.vector(gradient),
              Hessian = clm.hess(rho), ## ensure hessian is evaluated
              ## at optimum
              nll = nll, ## negative log-likelihood
              convergence = conv,
              ## 0: successful convergence
              ## 1: iteration limit reached
              ## 2: step factor reducd below minimum
              message = message,
              maxGradient = maxGrad,
              niter = c(outer = i, inner = innerIter),
              fitted = rho$fitted)
  return(res)
}


clm.nll <- function(rho) { ## negative log-likelihood
  with(rho, {
    eta1 <- drop(B1 %*% par) + o1
    eta2 <- drop(B2 %*% par) + o2
    fitted <- pfun(eta1) - pfun(eta2)
    if(all(fitted > 0))
### FIXME: rather if(all(is.finite(fitted)) && all(fitted > 0)) ???
      -sum(wts * log(fitted))
    else Inf
  })
}

clm.grad <- function(rho) { ## gradient of the negative log-likelihood
### return: vector of gradients - NA if fitted are not defined
  with(rho, {
    p1 <- dfun(eta1)
    p2 <- dfun(eta2)
    wtpr <- wts/fitted
    if(all(fitted > 0))
      -crossprod((B1 * p1 - B2 * p2), wtpr)
    else rep(NA, nalpha + nbeta)
  })
}

clm.hess <- function(rho) { ## hessian of the negative log-likelihood
### return Hessian matrix - NA if fitted are not defined
  with(rho, {
    dS.psi <- crossprod(B1 * gfun(eta1) * wtpr, B1) -
      crossprod(B2 * gfun(eta2) * wtpr, B2)
    dpi.psi <- B1 * p1 - B2 * p2
    if (all(fitted > 0))
      -dS.psi + crossprod(dpi.psi, (dpi.psi * wtpr / fitted))
    else array(NA, dim = c(nalpha + nbeta, nalpha + nbeta))
  })
}

Trace <- function(iter, stepFactor, val, maxGrad, par, first=FALSE) {
    t1 <- sprintf(" %3d:  %-5e:    %.3f:   %1.3e:  ",
                  iter, stepFactor, val, maxGrad)
    t2 <- formatC(par)
    if(first)
        cat("iter:  step factor:     Value:     max|grad|:   Parameters:\n")
    cat(t1, t2, "\n")
}

makeThresholds <- function(y, threshold) {
### Generate the threshold structure summarized in the transpose of
### the Jacobian matrix, tJac. Also generating nalpha and alphaNames.

### args:
### y - response variable, a factor
### threshold - one of "flexible", "symmetric" or "equidistant"

  stopifnot(is.factor(y))
  lev <- levels(y)
  ntheta <- nlevels(y) - 1
  
  if(threshold == "flexible") {
    tJac <- diag(ntheta)
    nalpha <- ntheta
    alphaNames <- paste(lev[-length(lev)], lev[-1], sep="|")
  }
  
  if(threshold == "symmetric") {
    if(!ntheta >=2)
      stop("symmetric thresholds are only meaningful for responses with 3 or more levels")
    if(ntheta %% 2) { ## ntheta is odd
      nalpha <- (ntheta + 1)/2 ## No. threshold parameters
      tJac <- t(cbind(diag(-1, nalpha)[nalpha:1, 1:(nalpha-1)],
                      diag(nalpha)))
      tJac[,1] <- 1
      alphaNames <-
        c("central", paste("spacing.", 1:(nalpha-1), sep=""))
    }
    else { ## ntheta is even
      nalpha <- (ntheta + 2)/2
      tJac <- cbind(rep(1:0, each = ntheta / 2),
                    rbind(diag(-1, ntheta / 2)[(ntheta / 2):1,],
                          diag(ntheta / 2)))
      tJac[,2] <- rep(0:1, each = ntheta / 2)
      alphaNames <- c("central.1", "central.2",
                      paste("spacing.", 1:(nalpha-2), sep=""))
    }
  }
  
  if(threshold == "equidistant") {
    if(!ntheta >=2)
      stop("equidistant thresholds are only meaningful for responses with 3 or more levels")
    tJac <- cbind(1, 0:(ntheta-1))
    nalpha <- 2
    alphaNames <- c("threshold.1", "spacing")
  }
  return(list(tJac = tJac, nalpha = nalpha, alphaNames = alphaNames))
}

clm.control <-
    function(method = c("Newton", "model.frame"), ## ...,
             trace = 0, maxIter = 100,
             gradTol = 1e-6, maxLineIter = 10)
{
### FIXME: change Newton to clm.fit?
  
  if(!all(is.numeric(c(maxIter, gradTol, maxLineIter))))
    stop("maxIter, gradTol, maxLineIter, convTol should all be numeric")
  
  ## method <- match.arg(method)
  ctrl <- list(method = match.arg(method),
               trace = as.integer(trace),
               maxIter = as.integer(maxIter),
               gradTol = as.numeric(gradTol),
               maxLineIter = as.integer(maxLineIter))
  
  ## if(convTol <= 0)
  ##   stop("convTol should be > 0")
  ## if(method == "Newton" && convTol > gradTol)
  ##   stop("convTol should be <= gradTol")
  
  ## list(method = method, ctrl = ctrl)
  return(ctrl)
}

