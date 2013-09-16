##' Coefficients (columns) are dropped from a design matrix to
##' ensure that it has full rank.
##'
##'  Redundant columns of the design matrix are identified with the
##'  LINPACK implementation of the \code{\link{qr}} decomposition and
##'  removed. The returned design matrix will have \code{qr(X)$rank}
##'  columns.
##' @title Ensure full rank design matrix
##' @param X     a design matrix, e.g., the result of
##' \code{\link{model.matrix}} possibly of less than full column rank,
##' i.e., with redundant parameters. Works for \code{ncol(X) >= 0} and
##' \code{nrow(X) >= 0}.
##' @param silent should a message not be issued if X is column rank
##' deficient?
##' @return The design matrix \code{X} without redundant columns.
##' @seealso \code{\link{qr}} and \code{\link{lm}}
##' @export
##' @keywords models
##' @author Rune Haubo Bojesen Christensen
##' @example inst/examples/examples-dropcoef.R
drop.coef <- function(X, silent = FALSE)
### works if ncol(X) >= 0 and nrow(X) >= 0
{
  ## test and match arguments:
  stopifnot(is.matrix(X))
  silent <- as.logical(silent)[1]
  ## perform the qr-decomposition of X using LINPACK methods:
  qr.X <- qr(X, tol = 1e-7, LAPACK = FALSE)
  if(qr.X$rank == ncol(X))
    return(X) ## return X if X has full column rank
  if(!silent) ## message the no. dropped columns:
    message(gettextf("design is column rank deficient so dropping %d coef",
                     ncol(X) - qr.X$rank))
  ## return the columns correponding to the first qr.x$rank pivot
  ## elements of X:
  newX <- X[, qr.X$pivot[1:qr.X$rank], drop = FALSE]
  ## did we succeed? stop-if-not:
  if(qr.X$rank != qr(newX)$rank)
    stop(gettextf("determination of full column rank design matrix failed"),
         call. = FALSE)
  return(newX)
}

