clm.control <-
  function(method = c("Newton", "model.frame", "ucminf", "nlminb",
             "optim"), ...,  trace = 0L, maxIter = 100L, gradTol = 1e-6,
           maxLineIter = 15L, relTol = 1e-6, tol = sqrt(.Machine$double.eps),
           maxModIter = 5L,
           convergence=c("warn", "silent", "stop", "message"))
{
  method <- match.arg(method)
  convergence <- match.arg(convergence)

  if(!all(is.numeric(c(maxIter, gradTol, maxLineIter, relTol, tol,
                       maxModIter))))
      stop("maxIter, gradTol, relTol, tol, maxModIter and maxLineIter should all be numeric")

  ctrl <- list(method = method,
               convergence = convergence,
               trace = as.integer(trace),
               maxIter = as.integer(maxIter),
               gradTol = as.numeric(gradTol),
               relTol = as.numeric(relTol),
               tol = as.numeric(tol),
               maxLineIter = as.integer(maxLineIter),
               maxModIter = as.integer(maxModIter))
if(method %in% c("ucminf", "nlminb", "optim"))
    ctrl$ctrl <- list(trace = as.integer(abs(trace)), ...)

  return(ctrl)
}

clmm.control <-
  function(method = c("ucminf", "nlminb", "model.frame"),
           ..., trace = 0, maxIter = 50, gradTol = 1e-4,
           maxLineIter = 50, useMatrix = FALSE,
           innerCtrl = c("warnOnly", "noWarn", "giveError"))
{
  method <- match.arg(method)
  innerCtrl <- match.arg(innerCtrl)
  useMatrix <- as.logical(useMatrix)
  stopifnot(is.logical(useMatrix))
  ctrl <- list(trace=if(trace < 0) 1 else 0,
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

  list(method = method, useMatrix = useMatrix, ctrl = ctrl,
       optCtrl = optCtrl)
}

## getCtrlArgs <- function(control, extras) {
## ### Recover control arguments from clmm.control and extras (...):
## ###
##   ## Collect control arguments in list:
##   ctrl.args <- c(extras, control$method, control$useMatrix,
##                  control$ctrl, control$optCtrl)
##   ## Identify the two occurences "trace", delete them, and add trace=1
##   ## or trace=-1 to the list of arguments:
##   which.trace <- which(names(ctrl.args) == "trace")
##   trace.sum <- sum(unlist(ctrl.args[which.trace]))
##   ctrl.args <- ctrl.args[-which.trace]
##   ## remove duplicated arguments:
##   ctrl.args <- ctrl.args[!duplicated(names(ctrl.args))]
##   if(trace.sum >= 1) ctrl.args$trace <- 1
##   if(trace.sum >= 2 || trace.sum <= -1) ctrl.args$trace <- -1
##   ## return the updated list of control parameters:
##   do.call("clmm.control", ctrl.args)
## }

getCtrlArgs <- function(control, extras) {
### Recover control arguments from clmm.control and extras (...):
###
    if(!is.list(control))
        stop("'control' should be a list")
    ## Collect control arguments in list:
    ## 1) assuming 'control' is a call to clmm.control:
        ctrl.args <-
        if(setequal(names(control), names(clmm.control())))
            c(extras, control["method"], control["useMatrix"],
              control$ctrl, control$optCtrl)
    ## assuming 'control' is specified with control=list( 'args'):
        else
            c(extras, control)
### NOTE: having c(extras, control) rather than c(control, extras)
### means that extras have precedence over control.
    ## Identify the two occurences "trace", delete them, and add trace=1
    ## or trace=-1 to the list of arguments:
    which.trace <- which(names(ctrl.args) == "trace")
    trace.sum <- sum(unlist(ctrl.args[which.trace]))
    if(trace.sum)
        ctrl.args <- ctrl.args[-which.trace]
    ## remove duplicated arguments:
    ctrl.args <- ctrl.args[!duplicated(names(ctrl.args))]
    if(trace.sum >= 1) ctrl.args$trace <- 1
    if(trace.sum >= 2 || trace.sum <= -1) ctrl.args$trace <- -1
    ## return the updated list of control parameters:
    do.call("clmm.control", ctrl.args)
}
