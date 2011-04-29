slice <- function(object, ...) {
  UseMethod("slice")
}

slice.clm <-
  function(object, which = seq_along(par), lambda = 3, grid = 1e2,
           quad.approx = TRUE, ...)
{
  ## argument matching and testing:
  stopifnot(is.numeric(lambda) && lambda > 0)
  stopifnot(is.numeric(grid) && grid >= 1)
  grid <- as.integer(grid)
  par <- coef(object)
  parNames <- names(par)
  npar <- length(par)
  stopifnot(length(which) == length(unique(which)))
  if(is.character(which))
    which <- match(which, parNames, nomatch = 0)
  if(!all(which %in% seq_along(par)))
    stop("invalid 'which' argument")
  stopifnot(length(which) > 0)
  which <- as.integer(which)
  ml <- object$logLik
  which.names <- parNames[which]

  ## get environment corresponding to object:
  rho <- update(object, doFit = FALSE)
  names(par) <- NULL
  rho$par <- par ## set rho$par to mle
  stopifnot(isTRUE(all.equal(clm.nll(rho), -object$logLik)))

  curv <- 1/diag(object$Hess) ## curvature in nll wrt. par
  par.range <- par + sqrt(curv) %o% c(-lambda, lambda)
  ## par.seq - list of length npar:
  par.seq <- sapply(which, function(ind) {
    seq(par.range[ind, 1], par.range[ind, 2], length = grid) }, 
                    simplify = FALSE)
  ## compute relative nll for all par.seq for each par:
  nll <- lapply(seq_along(which), function(i) { # for each par
    rho$par <- par ## reset par values to MLE
    sapply(par.seq[[ i ]], function(par.val) { # for each val
      rho$par[ which[i] ] <- par.val
      clm.nll(rho) + ml ## relative nll
    })
  })
  
  ## collect results in a list of data.frames:
  res <- lapply(seq_along(which), function(i) {
    structure(data.frame(par.seq[[ i ]], nll[[ i ]]),
              names = c(which.names[i], "nll"))
  })

  ## set attributes:
  names(res) <- which.names
  attr(res, "original.fit") <- object
  class(res) <- "slice.clm"

  if(!quad.approx) return(res)
  ## compute quadratic approx to nll:
  Quad <- function(par, mle, curv)
    ((mle - par)^2 / curv / 2)
  for(i in seq_along(which))
    res[[ i ]]$quad <-
      Quad(par.seq[[ i ]], par[ which[i] ], curv[ which[i] ])
  
  return(res)
}

plot.slice.clm <-
  function(x, which = seq_along(x), type = c("quadratic", "linear"),
           plot.mle = TRUE, ...)
{
  type <- match.arg(type)
  stopifnot(is.numeric(which))
  which <- as.integer(which)
  of <- attr(x, "original.fit")
  par <- coef(of)
  ml <- of$logLik

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
      abline(v = par[names(x)[i]])
  }
  
  return(invisible())
}
