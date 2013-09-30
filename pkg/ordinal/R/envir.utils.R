setLinks <- function(rho, link) {
### The Aranda-Ordaz and log-gamma links are not supported in this
### version of clm.
  rho$pfun <- switch(link,
                     logit = plogis,
                     probit = pnorm,
                     cloglog = function(x, lower.tail=TRUE) pgumbel(x,
                                             lower.tail=lower.tail, max=FALSE),
                     cauchit = pcauchy,
                     loglog = pgumbel,
                     "Aranda-Ordaz" = function(x, lambda) pAO(x, lambda),
                     "log-gamma" = function(x, lambda) plgamma(x, lambda))
  rho$dfun <- switch(link,
                     logit = dlogis,
                     probit = dnorm,
                     cloglog = function(x) dgumbel(x, max=FALSE),
                     cauchit = dcauchy,
                     loglog = dgumbel,
                     "Aranda-Ordaz" = function(x, lambda) dAO(x, lambda),
                     "log-gamma" = function(x, lambda) dlgamma(x, lambda))
  rho$gfun <- switch(link,
                     logit = glogis,
                     probit = gnorm,
                     cloglog = function(x) ggumbel(x, max=FALSE),
                     loglog = ggumbel,
                     cauchit = gcauchy,
                     "Aranda-Ordaz" = function(x, lambda) gAO(x, lambda), ## shouldn't happen
                     "log-gamma" = function(x, lambda) glgamma(x, lambda)
                     )
  rho$link <- link
}

makeThresholds <- function(ylevels, threshold) { ## , tJac) {
### Generate the threshold structure summarized in the transpose of
### the Jacobian matrix, tJac. Also generating nalpha and alpha.names.

### args:
### y - response variable, a factor
### threshold - one of "flexible", "symmetric" or "equidistant"
  ## stopifnot(is.factor(y))
  lev <- ylevels
  ntheta <- length(lev) - 1

  ## if(!is.null(tJac)) {
  ##   stopifnot(nrow(tJac) == ntheta)
  ##   nalpha <- ncol(tJac)
  ##   alpha.names <- colnames(tJac)
  ##   if(is.null(alpha.names) || anyDuplicated(alpha.names))
  ##     alpha.names <- as.character(1:nalpha)
  ##   dimnames(tJac) <- NULL
  ## }
  ## else { ## threshold structure identified by threshold argument:
    if(threshold == "flexible") {
      tJac <- diag(ntheta)
      nalpha <- ntheta
      alpha.names <- paste(lev[-length(lev)], lev[-1], sep="|")
    }

    if(threshold == "symmetric") {
      if(!ntheta >=2)
        stop("symmetric thresholds are only meaningful for responses with 3 or more levels",
             call.=FALSE)
      if(ntheta %% 2) { ## ntheta is odd
        nalpha <- (ntheta + 1)/2 ## No. threshold parameters
        tJac <- t(cbind(diag(-1, nalpha)[nalpha:1, 1:(nalpha-1)],
                        diag(nalpha)))
        tJac[,1] <- 1
        alpha.names <-
          c("central", paste("spacing.", 1:(nalpha-1), sep=""))
      }
      else { ## ntheta is even
        nalpha <- (ntheta + 2)/2
        tJac <- cbind(rep(1:0, each = ntheta / 2),
                      rbind(diag(-1, ntheta / 2)[(ntheta / 2):1,],
                            diag(ntheta / 2)))
        tJac[,2] <- rep(0:1, each = ntheta / 2)
        alpha.names <- c("central.1", "central.2",
                         paste("spacing.", 1:(nalpha-2), sep=""))
      }
    }
    ## Assumes latent mean is zero:
    if(threshold == "symmetric2") {
      if(!ntheta >=2)
        stop("symmetric thresholds are only meaningful for responses with 3 or more levels",
             call.=FALSE)
      if(ntheta %% 2) { ## ntheta is odd
        nalpha <- (ntheta - 1)/2 ## No. threshold parameters
        tJac <- rbind(apply(-diag(nalpha), 1, rev),
                      rep(0, nalpha),
                      diag(nalpha))
      }
      else { ## ntheta is even
        nalpha <- ntheta/2
        tJac <- rbind(apply(-diag(nalpha), 1, rev),
                      diag(nalpha))
      }
      alpha.names <- paste("spacing.", 1:nalpha, sep="")
    }

    if(threshold == "equidistant") {
      if(!ntheta >=2)
        stop("equidistant thresholds are only meaningful for responses with 3 or more levels",
             call.=FALSE)
      tJac <- cbind(1, 0:(ntheta-1))
      nalpha <- 2
      alpha.names <- c("threshold.1", "spacing")
    }
  ## }
  return(list(tJac = tJac, nalpha = nalpha, alpha.names = alpha.names))
}

