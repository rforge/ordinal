## This file contains:
## methods for computing, manipulating and extracting design matrices,
## weights, offsets, model.frames and things like that.

##  #################################
##  ## call sequence
##  clm() {
##      get_clmFormulas()
##      get_clm.mf()
##      get_clmTerms() # optionally
##      get_clmDesign()
##
##      makeThresholds()
##      drop.cols()
##
##      clm.newRho()
##      set.start()
##      setLinks()
##
##      clm.fit.NR()
##
##      clm.finalize()
##      extractFromFrames()
##      formatTheta()
##      conv.check()
##      get_clmInfoTab()
##  }
##
##  get_clmFormulas() {
##      getFullForm()
##  }
##
##  get_clm.mf() {
##      model.frame()
##  }
##
##  get_clmTerms() {
##      get_clm.mf()
##  }
##
##  get_clmDesign() {
##      checkContrasts()
##      get_clmDM() ## for formula, scale, nominal
##      getWeights()
##      get_clmY()
##  }
##
##  get_clmDM() {
##      model.matrix()
##      getContrasts()
##      getOffset()
##  }

get_clmTerms <-
    function(mc, formulas, call.envir=parent.frame(2L))
### Compute terms objects for each of the formulas.
{
    ## We need this approach in order to get the predvars and
    ## dataClasses attributes of the terms objects.
    nms <- c("formula", "scale", "nominal")
    keep <- match(nms, names(formulas), nomatch=0)
    lapply(formulas[keep], function(form) {
        terms(get_clm.mf(mc, form, attr(formulas, "envir"), call.envir))
    })
}

get_clmDesign <- function(fullmf, terms.list, contrasts) {
### Compute (y, X, wts, off, S, NOM etc.) for a clm object.
### clm-internal+external
###
### terms.list: list of terms.objects.
    stopifnot(all(sapply(terms.list, inherits, "terms")))

    ## Check that contrasts are specified only for terms in the model:
    checkContrasts(terms=attr(fullmf, "terms"), contrasts=contrasts)

    ## Extract X (design matrix for location effects) + terms, offset:
    res <- get_clmDM(fullmf, terms.list[["formula"]], contrasts,
                     type="formula")
    res$terms <- terms.list[["formula"]]
    res$off <- res$offset
    res$offset <- NULL

    ## Extract weights:
    res$wts <- getWeights(fullmf)

    ## Extract model response:
    res <- c(get_clmY(fullmf, res$wts), res)

    ## Extract S (design matrix for the scale effects):
    if(!is.null(terms.list$scale)) {
        ans <- get_clmDM(fullmf, terms.list[["scale"]], contrasts,
                         type="scale")
        res$S <- ans$X
        res$S.terms <- terms.list[["scale"]]
        res$S.off <- ans$offset
        if(attr(res$S.terms, "response") != 0)
            stop("response not allowed in 'scale'", call.=FALSE)
    }

    ## Extract NOM (design matrix for the nominal effects):
    if(!is.null(terms.list$nominal)) {
        ans <- get_clmDM(fullmf, terms.list[["nominal"]], contrasts,
                         type="nominal")
        res$NOM <- ans$X
        res$nom.terms <- terms.list[["nominal"]]
        if(attr(res$nom.terms, "response") != 0)
            stop("response not allowed in 'nominal'", call.=FALSE)
        if(!is.null(attr(res$nom.terms, "offset")))
            stop("offset not allowed in 'nominal'", call.=FALSE)
    }

    ## Return results (list of design matrices etc.):
    res
### NOTE: X, S and NOM are with dimnames and intercepts are
### guaranteed. They may be column rank deficient.
}

get_clmDM <-
    function(fullmf, terms, contrasts, check.intercept=TRUE,
             type="formula", get.offset=TRUE)
### Get DM (=Design Matrix):
{
    X <- model.matrix(terms, data=fullmf,
                      contrasts.arg=getContrasts(terms, contrasts))
    ## Test for intercept in X(?):
    Xint <- match("(Intercept)", colnames(X), nomatch = 0)
    if(check.intercept && Xint <= 0) {
        X <- cbind("(Intercept)" = rep(1, nrow(X)), X)
        warning(gettextf("an intercept is needed and assumed in '%s'", type),
                call.=FALSE)
    } ## Intercept in X is guaranteed.
    res <- list(X=X)
    if(get.offset)
        res$offset <- getOffset(fullmf, terms)
    res
}

get_clm.mf <-
    function(mc, formula, form.envir, call.envir=parent.frame(2L))
### clm-internal
### Extract the model.frame from formula
### mc - matched call containing: data, subset, weights, na.action
{
    ## Extract the full model.frame(mf):
    m <- match(c("data", "subset", "weights", "na.action"),
               names(mc), 0)
    mfcall <- mc[c(1, m)]
    mfcall$formula <- formula
    mfcall$drop.unused.levels <- TRUE
    mfcall[[1]] <- as.name("model.frame")
    if(is.null(mfcall$data)) mfcall$data <- form.envir
    eval(mfcall, envir=call.envir)
}

get_clmY <- function(fullmf, wts) {
    y <- model.response(fullmf, "any") ## any storage mode
    if(is.null(y)) stop("'formula' needs a response", call.=FALSE)
    if(!is.factor(y)) stop("response in 'formula' needs to be a factor", call.=FALSE)
    ## ylevels are the levels of y with positive weights
    ylevels <- levels(droplevels(y[wts > 0]))
    ## check that y has at least two levels:
    if(length(ylevels) == 1L)
        stop(gettextf("response has only 1 level ('%s'); expecting at least two levels",
                      ylevels), call.=FALSE)
    if(!length(ylevels))
        stop("response should be a factor with at least two levels")
    ## return:
    list(y=y, ylevels=ylevels)
}

get_clmFormulas <- function(mc, envir=parent.frame(2L))
### clm-internal
### Extracts and returns a list of formulas needed for further processing.
### mc: matched call
### envir: environment in which mc is to be evaluated
{
    ## Collect all variables in a full formula:
    ## evaluate the formulae in the enviroment in which clm was called
    ## (parent.frame(2)) to get them evaluated properly:
    forms <- list(eval(mc$formula, envir=envir))
    if(!is.null(mc$scale)) forms$scale <- eval(mc$scale, envir=envir)
    if(!is.null(mc$nominal)) forms$nominal <- eval(mc$nominal, envir=envir)
    ## get the environment of the formula. If this does not have an
    ## enviroment (it could be a character), then use the parent frame.
    form.envir <-
        if(!is.null(env <- environment(forms[[1]]))) env else envir
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
    forms$fullForm <- do.call("getFullForm", forms)
### FIXME: do we really need to set this name?
    names(forms)[1] <- "formula"
    ## set environment of 'fullForm' to the environment of 'formula':
    attr(forms, "envir") <- environment(forms$fullForm) <- form.envir
    ## return:
    forms
}

get_clmRho <- function(object, ...) {
  UseMethod("get_clmRho")
}


get_clmRho.default <-
    function(object, terms.list, contrasts, link, threshold,
             parent=parent.frame(), start=NULL, ...)
### .default method(?)
### object: model.frame (fullmf) with all variables present
### terms.list: list of terms.objects for each of the formulas in the
### clm object.
{
    ## Get design matrices etc:
    design <- get_clmDesign(fullmf=object,
                            terms.list=terms.list,
                            contrasts=contrasts)
    ## Get threshold information:
    design$ths <- makeThresholds(design$ylevels, threshold)
    ## Drop columns for aliased coefs:
    design <- drop.cols(design, drop.scale=FALSE, silent=TRUE)
    ## Set envir rho with variables: B1, B2, o1, o2, wts, fitted:
    rho <- with(design, {
        clm.newRho(parent.frame(), y=y, X=X, NOM=design$NOM, S=design$S,
                   weights=wts, offset=off, S.offset=design$S.off,
                   tJac=ths$tJac)
    })
    ## Set and check starting values for the parameters:
    start <- set.start(rho, start=start, get.start=is.null(start),
                       threshold=threshold, link=link, frames=design)
    rho$par <- as.vector(start) ## remove attributes
    ## Set pfun, dfun and gfun in rho:
    setLinks(rho, link)
    ## Return:
    rho
}

get_clmRho.clm <-
    function(object, parent=parent.frame(), ...) {
### Safely generate the model environment from a model object.
    o <- object
    get_clmRho.default(object=model.frame(o),
                       terms.list=o$terms.list,
                       contrasts=o$contrasts, start=c(o$start), link=o$link,
                       threshold=o$threshold, parent=parent, ...)
}

## get_mfcall <- function(mc, envir=parent.frame(2)) {
##     m <- match(c("data", "subset", "weights", "na.action"),
##                names(mc), 0)
##     mf <- mc[c(1, m)]
##     ## mf$formula <- fullForm
##     mf$drop.unused.levels <- TRUE
##     mf[[1]] <- as.name("model.frame")
##     ## if(is.null(mf$data)) mf$data <- form.envir
##     list(mfcall=mf, envir=parent.frame(2))
## }

