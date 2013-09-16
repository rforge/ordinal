##' Gradients of common densities
##'
##' Gradients of common density functions in their standard forms, i.e., with
##' zero location (mean) and unit scale. These are implemented in C++ for speed
##' and care is taken that the correct results are provided for the argument
##' being \code{NA}, \code{NaN}, \code{Inf}, \code{-Inf} or just extremely
##' small or large.
##'
##' The gradients are given by:
##' \itemize{
##'   \item{gnorm: If \eqn{f(x)} is the normal density with mean 0 and
##'     spread 1, then the gradient is \deqn{f'(x) = -x f(x)}
##'   }
##'   \item{glogis: If \eqn{f(x)} is the logistic density with mean 0 and
##'     scale 1, then the gradient is
##'     \deqn{f'(x) = 2 \exp(-x)^2 (1 + \exp(-x))^{-3} -
##'       \exp(-x)(1+\exp(-x))^{-2}}{f'(x) = 2 exp(-x)^2 (1 + exp(-x))^{-3} -
##'       exp(-x)(1+exp(-x))^{-2}}
##'   }
##'   \item{pcauchy: If
##'     \eqn{f(x) = [\pi(1 + x^2)^2]^{-1}}{f(x) =1 / [pi (1 + x^2)^2]}
##'     is the cauchy density with mean 0 and scale 1, then the gradient
##'     is
##'     \deqn{f'(x) = -2x [\pi(1 + x^2)^2]^{-1}}{f'(x) = -2x / [pi (1 +
##'       x^2)^2]}
##'   }
##' }
##'
##' @aliases gfun gnorm glogis gcauchy
##' @param x numeric vector of quantiles.
##' @param mean numeric scalar; mean.
##' @param sd numeric scalar; standard deviation.
##' @return a numeric vector of gradients.
##' @export
##' @useDynLib Rufus
##' @author Rune Haubo Bojesen Christensen
##' @family distributions
##' @seealso Gradients of densities are also implemented for the extreme value
##' distribtion (\code{\link[=dgumbel]{gumbel}}) and the the log-gamma
##' distribution (\code{\link[=lgamma]{log-gamma}}).
##' @keywords distribution
##' @example inst/examples/examples-gfun.R
gnorm <- function(x, mean = 0, sd = 1)
    .Call("gnorm", x, mean, sd, PACKAGE="Rufus")

##' @rdname gnorm
##' @export
glogis <- function(x)
    .Call("glogis", x, PACKAGE="Rufus")

##' @rdname gnorm
##' @export
gcauchy <- function(x)
    .Call("gcauchy", x, PACKAGE="Rufus")
