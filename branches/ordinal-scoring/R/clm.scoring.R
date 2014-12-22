clmfs <-
    function(formula, data, weights, start, subset, offset, na.action,
             contrasts, model = TRUE, trace=FALSE, convtol=1e-6,
             maxiter=10L)
{
    nll <- function() { ## Par, linpred
        nllvec <- vapply(ivec, function(i) {
            prob <- diff(c(0, gamma[[i]], 1))
            sum(Y[i, ] * log(prob))
        }, numeric(1))
        -sum(nllvec)
    }
    prepare <- function() {## X, Par, pfun, dfun
        if(p > 0)
            linpred <<- drop(X[, -1, drop=FALSE] %*% Par[-isq])
        eta <<- lapply(linpred, function(lp) Par[isq] - lp)
        gamma <<- lapply(eta, pfun)
        alpha <<- lapply(gamma, function(g) diff(c(0, g)))
        Lambda <- lapply(eta, function(ea) diag(dfun(ea), q))
        invisible()
    }
    score <- function() {
        grad <- numeric(1)
        for(i in ivec) {
            grad <- grad +
                t(Xtilde[[i]]) %*% Lambda[[i]] %*% dalpha.dgamma %*%
                    Ualpha(alpha[[i]], Y[i, ])
        }
        grad
    }
    expected <- function() {
        EFI <- array(0, dim=c(q+p, q+p))
        for(i in ivec) {
            A <- t(Xtilde[[i]]) %*% Lambda[[i]] %*% dalpha.dgamma
            EFI <- EFI + (A %*% efi(alpha[[i]], Y[i, ]) %*% t(A))
        }
        EFI
    }
    ## Initial argument matching and testing:
    mc <- match.call(expand.dots = FALSE)
    ## link <- match.arg(link)
    if(missing(formula)) stop("Model needs a formula")
    if(missing(contrasts)) contrasts <- NULL

    ## Compute: y, X, wts, off, mf:
    if (missing(data)) data <- environment(formula)
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "na.action",
                 "offset"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    Y <- model.response(mf, "any") ## any storage mode
    if(is.factor(Y)) stop("response should be a numeric matrix")
    if(is.data.frame(Y)) Y <- as.matrix(Y)
    stopifnot(all(is.numeric(Y)), is.matrix(Y), ncol(Y) > 1, all(Y >= 0))
    if(any(max(abs(Y - round(Y))) > 0))
        warning("non-integer response counts rounded and coerced to integer")
    Y[] <- as.integer(round(Y))
    ## design matrix:
    mt <- attr(mf, "terms")
    X <- if(!is.empty.model(mt)) {
        model.matrix(mt, mf, contrasts)
    } else cbind("(Intercept)" = rep(1, NROW(Y)))
    ## Test for intercept in X:
    Xint <- match("(Intercept)", colnames(X), nomatch = 0)
    if(Xint <= 0) {
        X <- cbind("(Intercept)" = rep(1, NROW(Y)), X)
        warning("an intercept is needed and assumed in 'formula'",
                call.=FALSE)
    } ## intercept in X is guaranteed.
    wts <- model.weights(mf)
    off <- model.offset(mf)
    ## Links:
    pfun <- pnorm
    dfun <- dnorm
    ## Setup:
    q <- ncol(Y) - 1; p <- ncol(X) - 1; n <- nrow(X); isq <- seq_len(q)
    ivec <- seq_len(n)
    ## Set starting values:
    nstart <- q + p
    if(missing(start)) {
        start <- c(cumsum(rep(0.1, q)), rep(0, p))
    } else {
        stopifnot(length(start) == nstart)
    }
    Par <- start
    ## More setup:
    dalpha.dgamma <-
        local({B <- diag(q); B[row(B)+1 == col(B)] <- -1; B})
    Xtilde <- lapply(ivec, function(i) {
        ## if(q == 1) -X[i, , drop=FALSE] else
        xi <- -X[i, -1, drop=FALSE]
        if(q > 1) xi <- xi[rep(1, q), ]
        cbind(diag(q), xi)
    })
    linpred <- rep(0, n)
    if(p > 0)
        linpred <- drop(X[, -1, drop=FALSE] %*% Par[-isq])
    eta <- lapply(linpred, function(lp) Par[isq] - lp)
    gamma <- lapply(eta, pfun)
    alpha <- lapply(gamma, function(g) diff(c(0, g)))
    Lambda <- lapply(eta, function(ea) diag(dfun(ea), q))
    if(!is.finite(nll())) stop("Likelihood not finite at starting values")
    ## Fisher scoring:
    for(i in seq_len(maxiter)) {
        if(trace) cat(i-1, "par:", paste(format(Par, digits=2)), "\n")
### NOTE: trace gives parameter estimates after iteration i. After
### convergence, S and H are evaluated at the converged parameters
### after iteration i.
        prepare()
        S <- score()
        H <- expected()
        Lt <- chol(H)
### TODO: step-half here? check that convtol is positive?
        y <- backsolve(Lt, S, transpose=TRUE)
        step <- c(backsolve(Lt, y))
        ## step <- solve(H, S)
        if(max(abs(step)) < convtol)
            break
        Par <- Par + step
    }
    ## Gather results:
    fitted <- do.call(rbind, lapply(gamma, function(g) diff(c(0, g, 1))))
    expected <- fitted * rowSums(Y)
    ll.vec <- vapply(ivec, function(i)
                     dmultinom(Y[i, ], prob=fitted[i, ], log=TRUE),
                     numeric(1))
### TODO: Should we name: kernel -> logLik and logLik -> fullLogLik?
### That would ensure compatability with clm objects.
### TODO: Should we reserve 'Hessian' for the observed Fisher
### information? Then we could name the expected Fisher information,
### simply 'Fisher'.
    res <- list(niter=i, coefficients=c(Par), Hessian=H, vcov=solve(H),
                last.step=step, nobs=prod(dim(Y)), logLik=sum(ll.vec),
                edf=p+q, xbeta=linpred, fitted=fitted, y=Y,
                expected=expected, kernel=-nll(), model=mf)
    other <- list(terms=mt, call=match.call(), contrasts=attr(X, "contrasts"),
                  xlevels=.getXlevels(mt, mf))
    ## Return:
    dev <- getDeviance(Y, kernel=res$kernel, edf=res$edf)
    res <- c(res, other, dev[c("null.deviance", "deviance", "df.null",
                               "df.residual")])
    class(res) <- "clmfs"
    res.names <- names(res)
    res[sort(res.names)]
}

clmfs2 <-
    function(formula, data, weights, start, subset, offset, na.action,
             contrasts, model = TRUE, trace=FALSE, convtol=1e-6,
             maxiter=10L)
{
    nll <- function(Y, X) {
        ## Compute log-likelihood for CLM with probit link
        if(p > 0)
            linpred <- drop(X[, -1, drop=FALSE] %*% Par[-isq])
        nllvec <- vapply(seq_along(linpred), function(i) {
            cumprob <- pfun(Par[isq] - linpred[i])
            prob <- diff(c(0, cumprob, 1))
            sum(Y[i, ] * log(prob))
        }, numeric(1))
        -sum(nllvec)
    }
    score <- function(Y, X) {
        if(p > 0)
            linpred <- drop(X[, -1, drop=FALSE] %*% Par[-isq])
        dalpha.dgamma <-
            local({B <- diag(q); B[row(B)+1 == col(B)] <- -1; B})
        grad <- lapply(seq_along(linpred), function(i) {
            eta <- Par[isq] - linpred[i]
            gamma <- pfun(eta) ## cumprob
            alpha <- diff(c(0, gamma))
            Lambda <- diag(dfun(eta))
            Xtilde <- cbind(diag(q), -X[i, -1])
            t(Xtilde) %*% Lambda %*% dalpha.dgamma %*% Ualpha(alpha, Y[i, ])
        })
        do.call("+", grad)
    }
    expected <- function(Y, X) {
        if(p > 0)
            linpred <- as.vector(X[, -1, drop=FALSE] %*% Par[-isq])
        dalpha.dgamma <-
            local({B <- diag(q); B[row(B)+1 == col(B)] <- -1; B})
        EFI <- lapply(seq_along(linpred), function(i) {
            eta <- Par[isq] - linpred[i]
            gamma <- pnorm(eta)
            alpha <- diff(c(0, gamma))
            Lambda <- diag(dfun(eta))
            Xtilde <- cbind(diag(q), -X[i, -1])
            A <- t(dalpha.dgamma) %*% Lambda %*% Xtilde
            t(A) %*% efi(alpha, Y[i, ]) %*% A
        })
        do.call("+", EFI)
    }
    ## Initial argument matching and testing:
    mc <- match.call(expand.dots = FALSE)
    ## link <- match.arg(link)
    if(missing(formula)) stop("Model needs a formula")
    if(missing(contrasts)) contrasts <- NULL

    ## Compute: Y, X, wts, off, mf:
    if (missing(data)) data <- environment(formula)
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "na.action",
                 "offset"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    Y <- model.response(mf, "any") ## any storage mode
    if(!is.matrix(Y)) Y <- as.matrix(Y)
    stopifnot(is.matrix(Y))
    ## design matrix:
    mt <- attr(mf, "terms")
    X <- if(!is.empty.model(mt)) {
        model.matrix(mt, mf, contrasts)
    } else cbind("(Intercept)" = rep(1, NROW(Y)))
    ## Test for intercept in X:
    Xint <- match("(Intercept)", colnames(X), nomatch = 0)
    if(Xint <= 0) {
        X <- cbind("(Intercept)" = rep(1, NROW(Y)), X)
        warning("an intercept is needed and assumed in 'formula'",
                call.=FALSE)
    } ## intercept in X is guaranteed.
    wts <- model.weights(mf)
    off <- model.offset(mf)
    ## Links:
    pfun <- pnorm
    dfun <- dnorm
    ## Setup:
    q <- ncol(Y) - 1
    p <- ncol(X) - 1
    isq <- seq_len(q)
    linpred <- rep(0, nrow(X))
    ## Set starting values:
    nstart <- q + p
    if(missing(start)) {
        start <- c(cumsum(rep(0.1, q)), rep(0, p))
    } else {
        stopifnot(length(start) == nstart)
    }
    Par <- start
    nll0 <- nll(Y, X)
    if(!is.finite(nll0)) stop("Likelihood not finite at starting values")
    ## Fisher scoring:
    for(i in seq_len(maxiter)) {
        if(trace) cat(i-1, "par:", paste(format(Par, digits=2)), "\n")
        S <- score(Y, X)
        H <- expected(Y, X)
        step <- solve(H, S)
        if(max(abs(step)) < convtol)
### FIXME: check that convtol is positive?
            break
        Par <- Par + step
    }
    res <- list(niter=i, coefficients=c(Par), Hessian=H, vcov=solve(H),
                last.step=step, nobs=prod(dim(Y)), logLik=nll(Y, X),
                edf=p+q)
    other <- list(terms=mt, call=match.call(), contrasts=attr(X, "contrasts"),
                  xlevels=.getXlevels(mt, mf))
    ## Return:
    c(res, other)
}

Ualpha <- function(alpha, n)  {
### Gradient/score function in alpha-notation.
    n <- as.vector(n)
    alpha <- as.vector(alpha)
    stopifnot(length(alpha) == length(n) - 1)
    J <- length(n)
    structure(n[-J] / alpha - n[J] / (1-sum(alpha)), dim=c(J-1, 1))
}
efi <- function(alpha, n) {
    J <- length(n)
    N <- sum(n)
    prob <- c(alpha, 1-sum(alpha))
    dpi.dalpha <- cbind(diag(J-1), -1)
    PI <- diag(sqrt(N / prob))
    tcrossprod(dpi.dalpha %*% PI)
    ## Alternatively:
    ## diag(sum(n)/alpha) + sum(n)/(1-sum(alpha))
}


#################################
## Manipulating factors to and from matrices in data.frames:
fac2mat <-
    function(formula, data=NULL, subset, na.action,
             stringsAsFactors=TRUE, table.name,
             zero.rows=c("keep", "remove"), ...)
{
    zero.rows <- match.arg(zero.rows)
    mc <- as.list(match.call(expand.dots=FALSE))
    mc$stringsAsFactors <- mc$table.name <- mc$... <- mc$zero.rows <- NULL
    ## Get ftable:
    mc[[1]] <- as.name("ftable")
    tab <- eval.parent(as.call(mc))
    colvar <- attr(tab, "col.vars")
    vars <- attr(tab, "row.vars")
    if(missing(table.name)) {
        table.name <- names(colvar)[1]
    } else if(table.name %in% names(vars)) {
        stop("'table.name' is already a variable name")
    }
    ## Make data.frame:
    df <- rev(expand.grid(rev(vars), stringsAsFactors=stringsAsFactors))
    mat <- as.matrix(tab)
    df[[table.name]] <- mat
    if(zero.rows == "remove") df <- df[rowSums(mat) > 0, ]
    df
}

mat2fac <-
    function(formula, data=NULL, subset, na.action,
             zero.freq=c("keep", "remove"), expand=FALSE, ...)
{
    zero.freq <- match.arg(zero.freq)
    mc <- match.call(expand.dots=FALSE)
    mc <- as.list(mc)
    mc$zero.freq <- mc$expand <- mc$... <- NULL
    mc[[1]] <- as.name("model.frame")
    mf <- eval.parent(as.call(mc))
    y <- model.response(mf)
    x <- mf[-attr(terms(mf), "response")]
    res <- mat2fac.mat(y, x)
    if(zero.freq == "remove") {
        res <- res[res$Freq > 0, ]
        rownames(res) <- seq_len(nrow(res))
    }
    if(expand) {
        res <- dfexpand(res, Freq="Freq")
    }
    res
}

mat2fac.mat <- function(mat, data) {
    stopifnot(is.matrix(mat))
    y <- structure(as.vector(mat), dim=dim(mat))
    stopifnot(all(sapply(data, is.factor)))
    stopifnot(prod(sapply(data, nlevels)) == nrow(mat))
    attr(y, "row.vars") <- lapply(data, levels)
    ## Attempt to extract
    dn <- dimnames(mat)
    if(!is.null(dn) && length(dn) == 2 && length(dn[[2]]) == ncol(mat)) {
        col.vars <- dn[2]
    } else {
        col.vars <- list(col.name = as.character(seq_len(ncol(mat))))
    }
    attr(y, "col.vars") <- col.vars
    class(y) <- "ftable"
    as.data.frame(y)
}

dfexpand <- function(data, Freq="Freq") {
    stopifnot(Freq %in% names(data))
    stopifnot(all(is.numeric(freq <- data[[Freq]])))
    data <- data[, -which(names(data) == Freq)]
    freq <- as.integer(round(freq))
    ind <- rep.int(seq_len(nrow(data)), freq)
    newdata <- data[ind, ]
    rownames(newdata) <- seq_len(nrow(newdata))
    newdata
}

dfreduce <- function(data, responseName="Freq", keep.all.levels=FALSE)
{
    stopifnot(is.data.frame(data))
    stopifnot(!(responseName %in% names(data)))
    uData <- unique(data)
    x <- apply(data, 1, paste, collapse="_")
    ux <- apply(uData, 1, paste, collapse="_")
    uData[[responseName]] <-
        as.vector(sapply(ux, function(u) sum(u == x)))
    uData <- uData[do.call(order, uData), ]
    rownames(uData) <- seq_len(nrow(uData))
    uData
}

## dfreduce2 <-
##     function(formula, data, subset, na.action, responseName="Freq")
## {
##     mc <- match.call(expand.dots=FALSE)
##     mc <- as.list(mc)
##     mc$responseName <- NULL
##     mc[[1]] <- as.name("model.frame")
##     mf <- eval.parent(as.call(mc))
##     ## Check that there is no right hand side.
##     ## Check that all left hand side variables are factors
##     ## Make expand.grid on all variables and their levels.
##     stopifnot(is.data.frame(data))
##     stopifnot(!(responseName %in% names(data)))
##     uData <- unique(data)
##     x <- apply(data, 1, paste, collapse="_")
##     uData[[responseName]] <- as.vector(table(x))
##     uData
## }
##

getDeviance <- function(Y, kernel, edf) {
### Compute deviance and degrees of freedom for null optionally
### current model.
###
### Y: response matrix
### kernel: (optional) log likelihood kernel for the current model
### edf: (optional) estimated degrees of freedom (number of
### parameters) in the current model.
    ## Checks:
    stopifnot(is.matrix(Y), all(is.numeric(Y)), all(Y >= 0))
    if(any(max(abs(Y - round(Y))) > 1e-8))
        warning("non-integer counts rounded and coerced to integer")
    Y[] <- as.integer(round(Y))
    ## Null/ log-likelihood:
    res <- vector("list")
    Z <- colSums(Y)
    null.pi <- Z / sum(Z)
    ## null log-likelihood kernel:
    res$null.kernel <- sum(Z * ifelse(null.pi > 0, log(null.pi), 0))
    ## log-likelihood kernel for saturated model:
    sat.pi <- Y / rowSums(Y)
    res$sat.kernel <- sum(Y * ifelse(sat.pi > 0, log(sat.pi), 0))
### Total (= null) deviance:
    res$null.deviance <- -2 * (res$null.kernel - res$sat.kernel)
    ## Degrees of Freedom:
    edf.sat <- nrow(Y) * (ncol(Y) - 1)
    edf.null <- ncol(Y) - 1
    res$df.null <- edf.sat - edf.null
    if(!missing(kernel))
        ## residual deviance is just the deviance:
        res$deviance <- -2 * (kernel - res$sat.kernel)
    if(!missing(edf))
        res$df.residual <- edf.sat - edf
    res
}
