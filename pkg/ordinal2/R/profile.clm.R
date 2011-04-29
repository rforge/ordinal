profile.clm2 <- ## using clm.fit()
  function(fitted, which.beta = seq_len(nbeta), alpha = 0.01,
           max.steps = 50, nsteps = 8, trace = FALSE,
           step.warn = 5, control = list(), ...)
### NOTE: seq_len(nbeta) works for nbeta = 0: numeric(0), while
### 1:nbeta gives c(1, 0).

### args:
### alpha - The likelihood is profiled in the 100*(1-alpha)%
###   confidence region as determined by the profile likelihood
### max.steps - the maximum number of profile steps in each direction
### nsteps - the approximate no. steps determined by the quadratic
### approximation to the log-likelihood function
### trace - if trace > 0 information of progress is printed
### step.warn - a warning is issued if the profile in each direction
###   contains less than step.warn steps (due to lack of precision).
{
  ## match and test arguments:
  stopifnot(is.numeric(alpha) && length(alpha) == 1 &&
            alpha > 0 && alpha < 1)
  stopifnot(round(max.steps) > round(nsteps))
  stopifnot(round(nsteps) > round(step.warn))
  stopifnot(round(nsteps) > 0 && round(step.warn) >= 0)
  max.steps <- round(max.steps)
  nsteps <- round(nsteps)
  step.warn <- round(step.warn)
  trace <- as.logical(trace)[1]
  beta.names <- names(fitted$beta)
  nbeta <- length(fitted$beta)
  if(is.character(which.beta)) 
    which.beta <- match(which.beta, beta.names, nomatch = 0)
  if(!all(which.beta %in% seq_len(nbeta)))
    stop("invalid 'parm' argument")
   stopifnot(length(which.beta) > 0)
  ## Extract various things from the original fit:
  orig.par <- coef(fitted) ## c(alpha, beta)
  beta0 <- fitted$beta ## regression coef.
  nalpha <- length(fitted$alpha) ## no. threshold coef.
  nbeta <- length(beta0)
  beta.names <- names(beta0)
  orig.logLik <- fitted$logLik
  std.err <- coef(summary(fitted))[-(1:nalpha), "Std. Error"]
  link <- fitted$link
  threshold <- fitted$threshold
  ## profile limit:
  lroot.max <- qnorm(1 - alpha/2)
  ## profile step length:
  delta <- lroot.max / nsteps
  ## results list:
  prof.list <- vector("list", length = length(which.beta)) 
  names(prof.list) <- beta.names[which.beta]
  ## get model.frame:
  mf <- update(fitted, method = "model.frame")
  y <- mf$y
  X <- mf$X
  wts <- mf$wts
  orig.off <- mf$off
  ## For each which.beta move up or down, fit the model and store the
  ## signed likelihood root statistic and parameter values:
  for(wb in which.beta) {
    par.wb <- matrix(orig.par, nrow = 1) ## MLE
    wb.name <- beta.names[wb]
    lroot.wb <- 0 ## lroot at MLE
    X.wb <- X[, -(1+wb), drop=FALSE] 
    for(direction in c(-1, 1)) { ## move down or up
      if(trace) {
        message("\nParameter: ", wb.name,
                c(" down", " up")[(direction + 1)/2 + 1])
        utils::flush.console()
      }
      ## (re)set starting values:
      start <- orig.par[-(nalpha + wb)]
      for(step in seq_len(max.steps)) {
        ## increment offset and refit model without wb parameter:
        beta.i <- beta0[wb] + direction * step * delta * std.err[wb]
        new.off <- orig.off + X[, 1+wb, drop=TRUE] * beta.i
        fit <- clm.fit(y=y, X=X.wb, weights=wts, offset=new.off,
                       control=control, start=start, link=link,
                       threshold=threshold)
        ## save likelihood root statistic:
        lroot <- -direction * sqrt(2*(orig.logLik + fit$nll))
        ## save lroot and pararameter values:
        lroot.wb <- c(lroot.wb, lroot)
        temp.par <- orig.par
        temp.par[names(fit$par)] <- fit$par
        temp.par[wb.name] <- beta.i
        par.wb <- rbind(par.wb, temp.par)
        ## update starting values:
        start <- fit$par
        ## break for loop if profile is far enough:
        if(abs(lroot) > lroot.max) break
      } ## end 'step in seq_len(max.steps)'
      ## test that lroot.max is reached and enough steps are taken: 
      if(abs(lroot) < lroot.max)
        warning("profile may be unreliable for ", wb.name,
                " because lroot.max was not reached for ",
                wb, c(" down", " up")[(direction + 1)/2 + 1])
      if(step <= step.warn)
        warning("profile may be unreliable for ", wb.name,
                " because only ", step, "\n  steps were taken ",
                c("down", "up")[(direction + 1)/2 + 1])
    } ## end 'direction in c(-1, 1)'
    ## order lroot and par. values and collect in a data.frame: 
    lroot.order <- order(lroot.wb, decreasing = TRUE)
    prof.list[[wb.name]] <-
      structure(data.frame(lroot.wb[lroot.order]), names = "lroot")
    prof.list[[wb.name]]$par.vals <- par.wb[lroot.order, ]
    
    if(!all(diff(par.wb[lroot.order, wb.name]) > 0))
      warning("likelihood is not monotonically decreasing from maximum,\n",
              "  so profile may be unreliable for ", wb.name)
  } ## end 'wb in which.beta'
  val <- structure(prof.list, original.fit = fitted) 
  class(val) <- c("profile.clm")
  return(val)
}

profile.clm <- ## using clm.fit.env()
  function(fitted, which.beta = seq_len(nbeta), alpha = 0.001,
           max.steps = 50, nsteps = 8, trace = FALSE,
           step.warn = 5, control = list(), ...)
### NOTE: seq_len(nbeta) works for nbeta = 0: numeric(0), while
### 1:nbeta gives c(1, 0).

### This is almost a copy of profile.clm2, which use clm.fit rather
### than clm.fit.env. The current implementation is the fastest, but
### possibly less readable.
{
  ## match and test arguments:
  stopifnot(is.numeric(alpha) && length(alpha) == 1 &&
            alpha > 0 && alpha < 1)
  stopifnot(round(max.steps) > round(nsteps))
  stopifnot(round(nsteps) > round(step.warn))
  stopifnot(round(nsteps) > 0 && round(step.warn) >= 0)
  max.steps <- round(max.steps)
  nsteps <- round(nsteps)
  step.warn <- round(step.warn)
  trace <- as.logical(trace)[1]
  beta.names <- names(fitted$beta)
  nbeta <- length(fitted$beta)
  if(is.character(which.beta)) 
    which.beta <- match(which.beta, beta.names, nomatch = 0)
  if(!all(which.beta %in% seq_len(nbeta)))
    stop("invalid 'parm' argument")
   stopifnot(length(which.beta) > 0)
  ## Extract various things from the original fit:
  orig.par <- coef(fitted) ## c(alpha, beta)
  beta0 <- fitted$beta ## regression coef.
  nalpha <- length(fitted$alpha) ## no. threshold coef.
  nbeta <- length(beta0)
  beta.names <- names(beta0)
  orig.logLik <- fitted$logLik
  std.err <- coef(summary(fitted))[-(1:nalpha), "Std. Error"]
  ## profile limit:
  lroot.max <- qnorm(1 - alpha/2)
  ## profile step length:
  delta <- lroot.max / nsteps
  ## results list:
  prof.list <- vector("list", length = length(which.beta)) 
  names(prof.list) <- beta.names[which.beta]
  ## get model.frame:
  mf <- update(fitted, method = "model.frame")
  X <- mf$X
  orig.off <- mf$off
  ## variables to set/update in the fitting environment, rho:
  B0 <- 1 * (col(matrix(0, length(mf$y), nlevels(mf$y))) ==
             c(unclass(mf$y)))
  o1 <- c(1e5 * B0[, ncol(B0)])
  o2 <- c(-1e5 * B0[, 1])
  rho <- update(fitted, doFit = FALSE)
  A1 <- rho$B1[, seq_len(nalpha), drop=FALSE]
  A2 <- rho$B2[, seq_len(nalpha), drop=FALSE]
  ## For each which.beta move up or down, fit the model and store the
  ## signed likelihood root statistic and parameter values:
  for(wb in which.beta) {
    par.wb <- matrix(orig.par, nrow = 1) ## MLE
    wb.name <- beta.names[wb]
    lroot.wb <- 0 ## lroot at MLE
    ## set variables in fitting environment:
    rho$B1 <- cbind(A1, -X[, -c(1, 1+wb), drop=FALSE])
    rho$B2 <- cbind(A2, -X[, -c(1, 1+wb), drop=FALSE])
    for(direction in c(-1, 1)) { ## move down or up
      if(trace) {
        message("\nParameter: ", wb.name,
                c(" down", " up")[(direction + 1)/2 + 1])
        utils::flush.console()
      }
      ## reset starting values:
      rho$par <- orig.par[-(nalpha+wb)]
      for(step in seq_len(max.steps)) {
        ## increment offset and refit model without wb parameter:
        beta.i <- beta0[wb] + direction * step * delta * std.err[wb]
        new.off <- orig.off + X[, 1+wb, drop=TRUE] * beta.i
        rho$o1 <- o1 - new.off
        rho$o2 <- o2 - new.off
        fit <- clm.fit.env(rho, control)
        ## save likelihood root statistic:
        lroot <- -direction * sqrt(2*(orig.logLik + fit$nll))
        ## save lroot and pararameter values:
        lroot.wb <- c(lroot.wb, lroot)
        temp.par <- orig.par
        temp.par[names(fit$par)] <- fit$par
        temp.par[wb.name] <- beta.i
        par.wb <- rbind(par.wb, temp.par)
        ## break for loop if profile is far enough:
        if(abs(lroot) > lroot.max) break
      } ## end 'step in seq_len(max.steps)'
      ## test that lroot.max is reached and enough steps are taken: 
      if(abs(lroot) < lroot.max)
        warning("profile may be unreliable for ", wb.name,
                " because lroot.max was not reached for ",
                wb, c(" down", " up")[(direction + 1)/2 + 1])
      if(step <= step.warn)
        warning("profile may be unreliable for ", wb.name,
                " because only ", step, "\n  steps were taken ",
                c("down", "up")[(direction + 1)/2 + 1])
    } ## end 'direction in c(-1, 1)'
    ## order lroot and par. values and collect in a data.frame: 
    lroot.order <- order(lroot.wb, decreasing = TRUE)
    prof.list[[wb.name]] <-
      structure(data.frame(lroot.wb[lroot.order]), names = "lroot")
    prof.list[[wb.name]]$par.vals <- par.wb[lroot.order, ]
    
    if(!all(diff(par.wb[lroot.order, wb.name]) > 0))
      warning("likelihood is not monotonically decreasing from maximum,\n",
              "  so profile may be unreliable for ", wb.name)
  } ## end 'wb in which.beta'
  val <- structure(prof.list, original.fit = fitted) 
  class(val) <- c("profile.clm")
  return(val)
}
  
confint.clm <-
  function(object, parm = seq_len(nbeta), level = 0.95,
           type = c("profile", "Wald"), trace = FALSE, ...)
{
  ## match and test arguments
  type <- match.arg(type)
  stopifnot(is.numeric(level) && length(level) == 1 &&
            level > 0 && level < 1)
  trace <- as.logical(trace)[1]
  nbeta <- length(object$beta)
  beta.names <- names(object$beta)
  if(is.character(parm)) 
    parm <- match(parm, beta.names, nomatch = 0)
  if(!all(parm %in% seq_len(nbeta)))
    stop("invalid 'parm' argument")
  stopifnot(length(parm) > 0)
  ## Wald CI:
  if(type == "Wald")
    return(confint.default(object = object, parm = beta.names[parm],
                           level = level))
  ## profile likelhood CI:
  message("Waiting for profiling to be done...")
  utils::flush.console()
  ## get profile:
  object <- profile(object, which.beta = beta.names[parm],
                    alpha = (1 - level)/4, trace = trace, ...)
  ## get and return CIs:
  confint(object, parm = beta.names[parm], level = level, ...)
}

confint.profile.clm <-
  function(object, parm = seq_len(nbeta), level = 0.95, ...)
{
  ## match and test arguments:
  stopifnot(is.numeric(level) && length(level) == 1 &&
            level > 0 && level < 1)
  of <- attr(object, "original.fit")
  beta.names <- names(object)
  nbeta <- length(beta.names)
  if(is.character(parm))
    parm <- match(parm, beta.names, nomatch = 0)
  if(!all(parm %in% seq_len(nbeta)))
    stop("invalid 'parm' argument")
  stopifnot(length(parm) > 0)
  ## prepare CI:
  a <- (1-level)/2
  a <- c(a, 1-a)
  pct <- paste(round(100*a, 1), "%")
  ci <- array(NA, dim = c(length(parm), 2),
              dimnames = list(beta.names[parm], pct))
  cutoff <- qnorm(a)
  ## compute CI from spline interpolation of the likelihood profile: 
  for(pm in parm) {
    pm <- beta.names[pm]
    pro <- object[[ pm ]]
    sp <- spline(x = pro[, "par.vals"][, pm], y = pro[, 1])
    ci[pm, ] <- approx(sp$y, sp$x, xout = rev(cutoff))$y
  }
  ## do not drop(ci) because rownames are lost for single coef cases:  
  return(ci) 
}

plot.profile.clm <-
    function(x, parm = seq_len(nbeta), level = c(0.95, 0.99),
             Log = FALSE, relative = TRUE, root = FALSE, fig = TRUE,
             approx = root, n = 1e3, ..., ylim = NULL)
{
  ## match and test arguments:
  stopifnot(is.numeric(level) && all(level > 0) &&
            all(level < 1))
  stopifnot(n == round(n) && n > 0)
  Log <- as.logical(Log)[1]
  relative <- as.logical(relative)[1]
  root <- as.logical(root)[1]
  fig <- as.logical(fig)[1]
  approx <- as.logical(approx)[1]
  of <- attr(x, "original.fit")
  beta <- of$beta
  beta.names <- names(x)
  nbeta <- length(beta.names)
  if(is.character(parm))
    parm <- match(parm, beta.names, nomatch = 0)
  if(!all(parm %in% seq_len(nbeta)))
    stop("invalid 'parm' argument")
  stopifnot(length(parm) > 0)
  ML <- of$logLik
  ## prepare return value:
  which.names <- beta.names[parm]
  spline.list <- vector("list", length(parm))
  names(spline.list) <- which.names
  if(approx) std.err <- sqrt(diag(vcov(of)))[which.names]
  ## for each pm make the appropriate plot:
  for(pm in beta.names[parm]) {
    ## confidence limits:
    lim <- sapply(level, function(x)
                  exp(-qchisq(x, df=1)/2) )
    pro <- x[[ pm ]]
    sp <- spline(x=pro[, "par.vals"][, pm], y=pro[, 1], n=n)
    if(approx) y.approx <- (beta[pm] - sp$x) / std.err[pm]
    if(root) {
      ylab <- "profile trace"
      lim <- c(-1, 1) %o% sqrt(-2 * log(lim))
      sp$y <- -sp$y
      if(approx) y.approx <- -y.approx
    } else { ## !root:
      sp$y <- -sp$y^2/2
      if(approx) y.approx <- -y.approx^2/2
      if(relative && !Log) {
        sp$y <- exp(sp$y)
        if(approx) y.approx <- exp(y.approx)
        ylab <- "Relative profile likelihood"
        if(missing(ylim)) ylim <- c(0, 1)
      }
      if(relative && Log) {
        ylab <- "Relative profile log-likelihood"
        lim <- log(lim)
      }
      if(!relative && Log) {
        sp$y <- sp$y + ML
        if(approx) y.approx <- y.approx + ML
        ylab <- "Profile log-likelihood"
        lim <- ML + log(lim)
      }
      if(!relative && !Log) 
        stop("Not supported: at least one of 'Log' and 'relative' ",
           "have to be TRUE")
    }
    spline.list[[ pm ]] <- sp
    if(fig) { ## do the plotting:
      plot(sp$x, sp$y, type = "l", ylim = ylim,
           xlab = pm, ylab = ylab, ...)
      abline(h = lim)
      if(approx) lines(sp$x, y.approx, lty = 2)
      if(root)  points(beta[pm], 0, pch = 3)
    }
  }
  attr(spline.list, "limits") <- lim
  invisible(spline.list)
}
