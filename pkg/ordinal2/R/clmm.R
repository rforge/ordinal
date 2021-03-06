clmm <-
  function(formula, data, weights, start, subset, 
           na.action, contrasts, Hess = TRUE, model = TRUE,
           link = c("logit", "probit", "cloglog", "loglog", 
             "cauchit"), ##, "Aranda-Ordaz", "log-gamma"), ## lambda, 
           doFit = TRUE, control,
           threshold = c("flexible", "symmetric", "equidistant"), ...)
{
### Extract the matched call and initial testing:
  mc <- match.call(expand.dots = FALSE)
  ## call clm() when there are no random effects.

  link <- match.arg(link)
  threshold <- match.arg(threshold)
  if(missing(formula))  stop("Model needs a formula")
  if(missing(contrasts)) contrasts <- NULL

  ## set control parameters:
  if(missing(control)) control <- clmm.control(...)
  if(!setequal(names(control), c("method", "ctrl", "optCtrl")))
    stop("specify 'control' via clmm.control()")

  ## Extract y, X, Zt, offset (off), weights (wts):
  frames <- clmm.model.frame(mc, contrasts)
  if(control$method == "model.frame") return(frames)

  ## Test rank deficiency and possibly drop some parameters:
  ## X is guarantied to have an intercept at this point.
  frames$X <- drop.coef(frames$X)

  ## Compute the transpose of the Jacobian for the threshold function,
  ## tJac and the names of the threshold parameters, alpha.names:
  ths <- makeThresholds(frames$y, threshold)

  ## set and update rho environment:
  rho <- with(frames, {
    clm.newRho(parent.frame(), y=y, X=X, weights=wts,
               offset=off, tJac=ths$tJac) })
  rho.clm2clmm(rho=rho, Zt=frames$Zt, grList=frames$grList,
               ctrl=control$ctrl)
  
### NOTE: alpha.names, beta.names, random.names, tau.names?
### nalpha, nbeta, nrandom, ntau
  nbeta <- rho$nbeta <- ncol(frames$X) - 1 ## no. fixef parameters
  nalpha <- rho$nalpha <- ths$nalpha ## no. threshold parameters
  ntau <- length(frames$grList) ## no. variance parameters
### FIXME: save ntau in rho?

  ## set starting values for the parameters:
  if(missing(start)) start <- clmm.start(frames, link, threshold)
  stopifnot(is.numeric(start) && 
            length(start) == (nalpha + nbeta + ntau))
  rho$par <- start
  
  ## Set pfun, dfun and gfun in rho:
  setLinks(rho, link)
  
  ## possibly return the environment, rho without fitting:
  if(!doFit)  return(rho)

  ## fit the clmm:
  fit <- clmm.fit.env(rho, control = control$optCtrl, Hess)
    
  ## Modify and return results:
  fit$link <- link
  fit$start <- start
  fit$threshold <- threshold
  fit$call <- match.call()
  fit$contrasts <- contrasts
  fit$tJac <- ths$tJac
  res <- clmm.finalize(fit=fit, frames=frames,
                       alpha.names=ths$alpha.names) 
  
  ## add model.frame to results list?
  if(model) res$model <- frames$mf
  
  return(res)
}

clmm.model.frame <- function(mc, contrasts) {
### mc - the matched call
### contrasts - contrasts for the fixed model terms

### FIXME: what if formula is not "evaluated", but just
### frm <- y ~ x; clmm2(frm, ...)?

  ## Construct a formula with all (fixed and random) variables
  ## (fullForm) and a formula with only fixed-effects variables
  ## (fixedForm):  
  fixedForm <- nobars(mc$formula) ## ignore terms with '|'
  fullForm <- subbars(mc$formula)      # substitute `+' for `|'
  
  ## Set the appropriate environments:
  environment(fullForm) <- environment(fixedForm) <-
    environment(eval(mc$formula))

  ## Extract full model.frame (fullmf):
  m <- match(c("data", "subset", "weights", "na.action", "offset"),
             names(mc), 0)
  mf <- mc[c(1, m)]
  mf$formula <- fullForm
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  fixedmf <- mf ## save call for later modification and evaluation
  fullmf <- eval(mf, envir = parent.frame(2)) ## '2' to get out of
  ## clmmFrames and clmm
### FIXME: What if data is a matrix?

  ## Extract model response:
  y <- model.response(fullmf)
  if(!is.factor(y)) stop("response needs to be a factor")
  
  ## Extract X:
  ## get X from fullmf to handle e.g., NAs and subset correctly
  fixedmf$formula <- fixedForm
  fixedmf <- eval(fixedmf, envir = parent.frame(2))
  fixedTerms <- attr(fixedmf, "terms")
  X <- model.matrix(fixedTerms, fullmf, contrasts)
  ## remove intercept from X:
  Xint <- match("(Intercept)", colnames(X), nomatch = 0)
  if(Xint <= 0) {
    X <- cbind("(Intercept)" = rep(1, n), X)
    warning("an intercept is needed and assumed")
  } ## intercept in X is garanteed.

  ## Make grList:
  barList <- expandSlash(findbars(mc$formula[[3]]))
  if (!length(barList))
    stop("No random effects terms specified in formula")
  names(barList) <- unlist(lapply(barList,
                                  function(x) deparse(x[[3]])))
  ## get list of grouping factors from model.frame:
  grList <- lapply(barList, function(x)
             eval(substitute(as.factor(fac)[,drop = TRUE],
                             list(fac = x[[3]])), fullmf))
  ## order random terms with decreasing no. levels:
  grList <- grList[rev(order(sapply(grList, nlevels)))]

  ## test that only random effects on the intercept are specified:
  grChar <- unlist(lapply(barList, function(x) as.character(x[[2]])))
  if(any(grChar != "1"))
    stop(gettextf("random terms has to be on the form '(1|factor)'"),
         call. = FALSE)
  
  ## test that all variables for the random effects are factors and
  ## have at least 3 levels: 
  stopifnot(all(sapply(grList, is.factor)))
  stopifnot(all(sapply(grList, nlevels) > 2))

  ## generate transpose of ranef model matrix:
  Zt <- do.call(rBind, lapply(grList, as, "sparseMatrix"))
  
  ## extract weights and offset:
  n <- nrow(X)
  if(is.null(wts <- model.weights(fullmf))) wts <- rep(1, n)
  if(is.null(off <- model.offset(fullmf))) off <- rep(0, n)
  
  ## check weights and offset:
  if (any(wts <= 0))
    stop(gettextf("negative weights or weights of zero are not allowed"))
  if(length(wts) && length(wts) != NROW(y))
    stop(gettextf("number of weights is %d should equal %d
(number of observations)", length(wts), NROW(y)))
  if(length(off) && length(off) != NROW(y))
    stop(gettextf("number of offsets is %d should equal %d
(number of observations)", length(off), NROW(y)))
  
  ## set fixef terms attribute of fullmf:
  attr(fullmf, "terms") <- fixedTerms
  res <- list(y = y, X = X, Zt = Zt, wts = as.double(wts),
              off = as.double(off), mf = fullmf, grList = grList)
  ## Note: X is with dimnames and an intercept is guarantied.
  return(res)
}

rho.clm2clmm <- function(rho, Zt, grList, ctrl)
### update environment, rho returned by clm.
{
### FIXME: write default list of control arguments?
  rho$ctrl = ctrl
  rho$nlev <- as.vector(sapply(grList, nlevels))
  rho$Zt <- Zt
  rho$random.names <- rownames(Zt)
  dimnames(rho$Zt) <- list(NULL, NULL)
  rho$nrandom <- sum(rho$nlev) ## no. random effects
  rho$ntau <- length(rho$nlev) ## no. random terms
  rho$tau <- rep(0, rho$ntau) ## with(rho, exp(par[nalpha + nbeta + 1:s]))
  rho$varVec <- rep.int(rho$tau, rho$nlev)
  rho$Vt <- crossprod(Diagonal(x = rho$varVec), rho$Zt)
  ## rho$Lambda <- Diagonal(x = rep.int(rho$tau, rho$nlev))
  ## rho$Vt <- crossprod(rho$Lambda, rho$Zt)
  rho$L <- Cholesky(tcrossprod(rho$Vt), LDL = TRUE, super = FALSE,
                    Imult = 1)
  rho$Niter <- 0
  rho$u <- rho$uStart <- rep(0, rho$nrandom)
  rho$.f <- if(package_version(packageDescription("Matrix")$Version) >
               "0.999375-30") 2 else 1
}

clmm.fit.env <-
  function(rho, control = list(), Hess = FALSE) {
  ## fit the model with Laplace:
  fit <- ucminf(rho$par, function(par) getNLA(rho, par),
                control = control)
  
  ## ensure evaluation at optimum:
  nllFast.u(rho)
  update.u(rho)

  ## prepare list of results:
  res <- list(coefficients = fit$par, optRes = fit,
              logLik = -fit$value, fitted = rho$fitted,
              random = rho$u, L = rho$L, Niter = rho$Niter)
### FIXME: set res$condVar

  ## compute hessian?
  if(Hess)
    res$Hessian <-
      hessian(function(par) getNLA(rho, par), x = fit$par,
              method.args = list(r = 2, show.details = TRUE))

  return(res)
}

getNLA <- function(rho, par) {
### negative log-likelihood by the Laplace approximation
    if(!missing(par)) rho$par <- par
    if(!update.u(rho)) return(Inf)
    if(any(rho$D < 0)) return(Inf)
    logDetD <- c(suppressWarnings(determinant(rho$L)$modulus)) -
      rho$nrandom * log(2*pi) / 2
    rho$nll + logDetD 
}

nll.u <- function(rho) { ## negative log-likelihood
### Not allowing for scale, and flexible link functions
  rho$tau <- with(rho, exp(par[nalpha + nbeta + 1:ntau]))
  rho$varVec <- rep.int(rho$tau, rho$nlev)
  ## rho$b <- rho$varVec * rho$u
  b.exploded <- as.vector(crossprod(rho$Zt, rho$varVec * rho$u))
  rho$eta1Fix <- drop(rho$B1 %*% rho$par[1:(rho$nalpha + rho$nbeta)])
  rho$eta2Fix <- drop(rho$B2 %*% rho$par[1:(rho$nalpha + rho$nbeta)])
### FIXME: possibly it should be '- b.exploded' ?
  rho$eta1 <- as.vector(rho$eta1Fix + b.exploded + rho$o1)
  rho$eta2 <- as.vector(rho$eta2Fix + b.exploded + rho$o2)
  rho$fitted <- rho$pfun(rho$eta1) - rho$pfun(rho$eta2)
  if(any(is.na(rho$fitted)) || any(rho$fitted <= 0))
    nll <- Inf
  else
    nll <- -sum(rho$wts * log(rho$fitted)) -
      sum(dnorm(x=rho$u, mean=0, sd=1, log=TRUE))
  nll
}


nllFast.u <- function(rho) { ## negative log-likelihood
### Not allowing for scale, and flexible link functions
  rho$varVec <- rep.int(rho$tau, rho$nlev)
  b.exploded <- as.vector(crossprod(rho$Zt, rho$varVec * rho$u))
### FIXME: possibly it should be '- b.exploded' ?
  rho$eta1 <- as.vector(rho$eta1Fix + b.exploded + rho$o1)
  rho$eta2 <- as.vector(rho$eta2Fix + b.exploded + rho$o2)
  rho$fitted <- rho$pfun(rho$eta1) - rho$pfun(rho$eta2)
  if(any(is.na(rho$fitted)) || any(rho$fitted <= 0))
    nll <- Inf
  else
    nll <- -sum(rho$wts * log(rho$fitted)) -
      sum(dnorm(x=rho$u, mean=0, sd=1, log=TRUE))
  nll
}


grad.u <- function(rho){
### should only be called with up to date values of eta1, eta2, par
  ## compute phi1:
  rho$p1 <- rho$dfun(rho$eta1)
  rho$p2 <- rho$dfun(rho$eta2)
  rho$wtpr <- rho$wts/rho$fitted 
  phi1 <- as.vector(rho$wtpr * (rho$p2 - rho$p1))
  return( (rho$Zt %*% phi1) * rho$varVec + rho$u )
}

hess.u <- function(rho) {
### should only be called with up-to-date values of eta1, eta2, par,
### p1, p2
  ## browser()
  g1 <- rho$gfun(rho$eta1) ## does not need to be saved
  g2 <- rho$gfun(rho$eta2) ## does not need to be saved
  phi2 <- ((rho$p1 - rho$p2)^2 / rho$fitted - g1 + g2) * rho$wtpr
  ## rho$Lambda <- Diagonal(x = rho$varVec)
  ## rho$Vt <- crossprod(rho$Lambda,
  ##                     tcrossprod(rho$Zt, Diagonal(x = sqrt(phi2))))
  rho$Vt <- crossprod(Diagonal(x = rho$varVec),
                      tcrossprod(rho$Zt, Diagonal(x = sqrt(phi2))))
  return(update(rho$L, rho$Vt, mult = 1))
  ## return(tcrossprod(rho$Vt) + Diagonal(rho$nrandom)) ## Hessian
}

update.u <- function(rho)
{
  stepFactor <- 1
  innerIter <- 0
  rho$u <- rho$uStart
  rho$nll <- nll.u(rho)
  if(!is.finite(rho$nll)) return(FALSE)
  rho$gradient <- grad.u(rho)
  maxGrad <- max(abs(rho$gradient))
  conv <- -1  ## Convergence flag
  message <- "iteration limit reached when updating the random effects"
  if(rho$ctrl$trace > 0)
    Trace(iter=0, stepFactor, rho$nll, maxGrad, rho$u, first=TRUE)
  ## Newton-Raphson algorithm:
  for(i in 1:rho$ctrl$maxIter) {
    if(maxGrad < rho$ctrl$gradTol) {
      message <- "max|gradient| < tol, so current
iterate is probably solution"
      if(rho$ctrl$trace > 0)
        cat("\nOptimizer converged! ", "max|grad|:",
            maxGrad, message, fill = TRUE)
            conv <- 0
      break
    }
    rho$L <- hess.u(rho)
    step <- as.vector(solve(rho$L, rho$gradient))
    rho$u <- rho$u - stepFactor * step
    nllTry <- nllFast.u(rho) ## no 'X %*% beta' update
    lineIter <- 0
    
    ## Step halfing:
    while(nllTry > rho$nll) {
      stepFactor <- stepFactor/2
      rho$u <- rho$u + stepFactor * step
      nllTry <- nllFast.u(rho) ## no 'X %*% beta' update
      lineIter <- lineIter + 1
      if(rho$ctrl$trace > 0)
        Trace(i+innerIter, stepFactor, rho$nll, maxGrad,
              rho$u, first=FALSE)
      if(lineIter > rho$ctrl$maxLineIter){
        message <- "step factor reduced below minimum when updating
the random effects"
        conv <- 1
        break
      }
      innerIter <- innerIter + 1
    }
    rho$nll <- nllTry
    rho$gradient <- grad.u(rho)
    maxGrad <- max(abs(rho$gradient))
    if(rho$ctrl$trace > 0)
      Trace(i+innerIter, stepFactor, rho$nll, maxGrad, rho$u, first=FALSE)
    stepFactor <- min(1, 2 * stepFactor)
  }
  if(conv != 0 && rho$ctrl$innerCtrl == "warnOnly") {
    warning(message, "\n  at iteration ", rho$Niter)
    utils::flush.console()
  }
  else if(conv != 0 && rho$ctrl$innerCtrl == "giveError")
        stop(message, "\n  at iteration ", rho$Niter)
  rho$Niter <- rho$Niter + i
    rho$L <- hess.u(rho)
  if(!is.finite(rho$nll))
    return(FALSE)
  else
    return(TRUE)
}

clmm.finalize <-
  function(fit, frames, alpha.names = ths$alpha.names)
{
  ## get coefficient names and lengths:
  beta.names <- colnames(frames$X)[-1]
  random.names <- sapply(frames$grList, levels)
  tau.names <- names(frames$grList)  
  nalpha <- length(alpha.names)
  nbeta <- length(beta.names)
  ntau <- length(tau.names)
  fit$nlev <- sapply(frames$grList, nlevels)

  ## test appropriate length of coefficients:
  stopifnot(length(fit$coefficients) == nalpha + nbeta + ntau)
  
  fit <- within(fit, {
    ## extract coefficients from 'fit':
    alpha <- coefficients[1:nalpha]
    names(alpha) <- alpha.names
    beta <- if(nbeta > 0) coefficients[nalpha + 1:nbeta]
    else numeric(0)
    names(beta) <- beta.names
    tau <- coefficients[nalpha + nbeta + 1:ntau]
    stDev <- exp(tau)
    names(stDev) <- tau.names
    names(tau) <- paste("log", tau.names, sep = ".")
    coefficients <- c(alpha, beta, tau)
    if(exists("Hessian", inherits = FALSE)) {
      dimnames(Hessian) <- list(names(coefficients),
                                names(coefficients))
    }
    varMat <- matrix(c(stDev^2, stDev),
                     nrow = length(stDev), ncol=2)
    dimnames(varMat) <- list(tau.names, c("Var", "Std.Dev"))
    ## set various fit elements:
    edf <- length(coefficients) ## estimated degrees of freedom
    nobs <- sum(frames$wts)
    n <- length(frames$wts)
    fitted.values <- fitted
    df.residual = nobs - edf
    terms <- frames$terms
    na.action <- frames$na.action
    info <-
      data.frame("link" = link,
                 "threshold" = threshold,
                 "nobs" = nobs, 
                 "logLik" = formatC(logLik, digits=2, format="f"),
                 "AIC" = formatC(-2*logLik + 2*edf, digits=2,
                   format="f"),
                 "niter" = paste(optRes$info["neval"], "(", Niter, ")",
                   sep=""), 
                 "max.grad" = formatC(optRes$info["maxgradient"], digits=2,
                   format="e") 
                 ## BIC is not part of output since it is not clear what
                 ## the no. observations are. 
                 )
    ## compute ranef estimates and conditional variance:
    ranef <- rep.int(stDev, nlev) * random
    names(ranef) <- as.vector(unlist(random.names))
    ranef <- split(x=-ranef, f=rep.int(names(frames$grList), nlev))
### NOTE: parameterization of random effects should change sign; this
### would avoid changing the sign at this point.
    condVar <- as.vector(diag(solve(L)) * rep.int(stDev^2, nlev))
    names(condVar) <- as.vector(unlist(random.names))
    condVar <- split(x=condVar, f=rep.int(names(frames$grList), nlev)) 

    ## remove excess elements:
    rm(list = c("L")) 
  })
  ## set class and return fit:
  class(fit) <- "clmm"
  return(fit)
}

clmm.control <-
  function(method = c("ucminf", "model.frame"),
           ..., trace = 0, maxIter = 50, gradTol = 1e-3,
           maxLineIter = 50,
           innerCtrl = c("warnOnly", "noWarn", "giveError"))
{
  method <- match.arg(method)
  innerCtrl <- match.arg(innerCtrl)
  ctrl <- list(trace=ifelse(trace < 0, 1, 0),
               maxIter=maxIter,
               gradTol=gradTol,
               maxLineIter=maxLineIter,
               innerCtrl=innerCtrl)
  optCtrl <- list(trace = abs(trace), ...)
  
  if(!is.numeric(unlist(ctrl[-5])))
    stop("maxIter, gradTol, maxLineIter and trace should all be numeric")
  if(any(ctrl[-c(1, 5)] <= 0))
    stop("maxIter, gradTol and maxLineIter have to be > 0")
  if(method == "ucminf" && !"grtol" %in% names(optCtrl))
    optCtrl$grtol <- 1e-5
  if(method == "ucminf" && !"grad" %in% names(optCtrl))
    optCtrl$grad <- "central"
  
  list(method = method, ctrl = ctrl, optCtrl = optCtrl)
}
