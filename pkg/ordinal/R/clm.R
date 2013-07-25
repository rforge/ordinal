clm.newRho <-
  function(parent, y, X, weights, offset, tJac)
### Set variables in rho: B1, B2, o1, o2 and wts.
{
  rho <- new.env(parent = parent)

  ## Make B1, B2, o1, o2 based on y, X and tJac:
  ntheta <- nlevels(y) - 1
  n <- nrow(X)
  B2 <- 1 * (col(matrix(0, n, ntheta + 1)) == c(unclass(y)))
  rho$o1 <- c(1e5 * B2[, ntheta + 1]) - offset
  rho$o2 <- c(-1e5 * B2[,1]) - offset
  B1 <- B2[, -(ntheta + 1), drop = FALSE]
  B2 <- B2[, -1, drop = FALSE]
  ## adjust B1 and B2 for structured thresholds:
  rho$B1 <- B1 %*% tJac
  rho$B2 <- B2 %*% tJac
  ## update B1 and B2 with location effects (X):
  nbeta <- NCOL(X) - 1
  if(nbeta > 0) {
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
  function(y, X, S, N, weights = rep(1, nrow(X)),
           offset = rep(0, nrow(X)), S.offset = rep(0, nrow(X)),
           control = list(), start,
           link = c("logit", "probit", "cloglog", "loglog", "cauchit"),
           threshold = c("flexible", "symmetric", "symmetric2", "equidistant"))
### This function basically does the same as clm, but without setting
### up the model matrices from formulae, and with minimal post
### processing after parameter estimation.
{
  ## Initial argument matching and testing:
  threshold <- match.arg(threshold)
  link <- match.arg(link)
  control <- do.call(clm.control, control)
  if(missing(y)) stop("please specify y")
  if(missing(X)) X <- cbind("(Intercept)" = rep(1, length(y)))
  stopifnot(is.factor(y), is.matrix(X))
  if(missing(weights) || is.null(weights))
      weights <- rep(1, length(y))
  if(missing(offset) || is.null(offset))
      offset <- rep(0, length(y))
  if(missing(S.offset) || is.null(S.offset))
      S.offset <- rep(0, length(y))
  stopifnot(length(y) == nrow(X) &&
            length(y) == length(weights) &&
            length(y) == length(offset) &&
            length(y) == length(S.offset))
  frames <- list(y=y, X=X)
  y[weights <= 0] <- NA
  frames$ylevels <- levels(droplevels(y))
  ## S and N are optional:
  if(!missing(S) && !is.null(S)) {
    frames$S <- S
    stopifnot(is.matrix(S),
              length(y) == nrow(S))
  }
  if(!missing(N) && !is.null(N)) {
    frames$NOM <- N
    stopifnot(is.matrix(N),
              length(y) == nrow(N))
  }

  ## Identify model as 'simple clm' or 'extended clm':
  Class <- if(any(S.offset != 0) || !missing(S) || link == "cauchit")
    c("eclm", "clm") else c("sclm", "clm")
  ## get  threshold structure:
  frames$ths <- makeThresholds(frames$ylevels, threshold)
  ## test for column rank deficiency in design matrices:
  frames <- drop.cols(frames, silent=TRUE)

  ## Set envir rho with variables: B1, B2, o1, o2, wts, fitted...:
  rho <- eclm.newRho(parent.frame(), y=frames$y, X=frames$X,
                     NOM=frames$NOM, S=frames$S, weights=weights,
                     offset=offset, S.offset=S.offset,
                     tJac=frames$ths$tJac)

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

  ## Set inverse link function and its derivatives (pfun, dfun and
  ## gfun) in rho:
  setLinks(rho, link)

  ## Fit the model:
  fit <- if(control$method == "Newton") {
      clm.fit.NR(rho, control) } else {
          clm.fit.optim(rho, control$method, control$ctrl) }

  ## Format and return the fit:
  fit$coef.names <- frames$coef.names
  fit$aliased <- lapply(frames$aliased, as.logical)
  if(control$method == "Newton" &&
     !is.null(start.iter <- attr(start, "start.iter")))
    fit$niter <- fit$niter + start.iter

  return(fit)
### FIXME: should 'par' be 'coefficients' to allow coef(fit) etc.?
### FIXME: should fit contain 'vcov'?
}

clm.nll <- function(rho) { ## negative log-likelihood
### For linear models
  with(rho, {
    eta1 <- drop(B1 %*% par) + o1
    eta2 <- drop(B2 %*% par) + o2
  })
### NOTE: getFitted is not found from within rho, so we have to
### evalueate it outside of rho
  rho$fitted <- getFittedC(rho$eta1, rho$eta2, rho$link)
  if(all(rho$fitted > 0))
### NOTE: Need test here because some fitted <= 0 if thresholds are
### not ordered increasingly.
### It is assumed that 'all(is.finite(pr)) == TRUE'
    -sum(rho$wts * log(rho$fitted))
  else Inf
}

clm.grad <- function(rho) { ## gradient of the negative log-likelihood
### return: vector of gradients
### For linear models
  with(rho, {
    p1 <- dfun(eta1)
    p2 <- dfun(eta2)
    wtpr <- wts/fitted
    dpi.psi <- B1 * p1 - B2 * p2
    -crossprod(dpi.psi, wtpr)
### NOTE: It is assumed that all(fitted > 0) == TRUE and that
### all(is.finite(c(p1, p2))) == TRUE
  })
}

clm.hess <- function(rho) { ## hessian of the negative log-likelihood
### return Hessian matrix
### For linear models
  with(rho, {
    dg.psi <- crossprod(B1 * gfun(eta1) * wtpr, B1) -
      crossprod(B2 * gfun(eta2) * wtpr, B2)
    -dg.psi + crossprod(dpi.psi, (dpi.psi * wtpr / fitted))
### NOTE: It is assumed that all(fitted > 0) == TRUE and that
### all(is.finite(c(g1, g2))) == TRUE
  })
}


