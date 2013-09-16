##' The Gumbel distribution
##'
##' Density, distribution function, quantile function, random generation, and
##' gradient of density of the extreme value (maximum and minimum)
##' distributions. The Gumbel distribution is also known as the extreme value
##' maximum distribution, the double-exponential distribution and the
##' log-Weibull distribution.
##'
##' \code{dgumbel}, \code{pgumbel} and \code{ggumbel} are implemented in C++ for
##' speed and care is taken that 'correct' results are provided for values of
##' \code{NA}, \code{NaN}, \code{Inf}, \code{-Inf} or just extremely small or
##' large.
##'
##' @aliases dgumbel pgumbel qgumbel rgumbel ggumbel, gumbel
##' @param q,x numeric vector of quantiles.
##' @param location numeric scalar.
##' @param scale numeric scalar.
##' @param lower.tail logical; if \code{TRUE} (default), probabilities are
##' \eqn{P[X \leq x]}{P[X <= x]} otherwise, \eqn{P[X > x]}.
##' @param max distribution for extreme maxima (default) or minima? The default
##' corresponds to the standard right-skew Gumbel distribution.
##' @return \code{pgumbel} gives the distribution function, \code{dgumbel}
##' gives the density, \code{ggumbel} gives the gradient of the density,
##' \code{qgumbel} is the quantile function, and \code{rgumbel} generates
##' random deviates.
##' @author Rune Haubo Bojesen Christensen
##' @seealso Gradients of densities are also implemented for the normal,
##' logistic, cauchy, cf. \code{\link[=gnorm]{gfun}} and the log-gamma
##' distribution, cf. \code{\link{lgamma}}.
##' @references \url{wikipedia.org/wiki/Gumbel_distribution}
##' @export
##' @useDynLib Rufus
##' @family distributions
##' @keywords distribution
##' @example inst/examples/examples-gumbel.R
pgumbel <-
    function(q, location = 0, scale = 1, lower.tail = TRUE, max = TRUE)
### Currently only unit length location and scale are supported.
    .Call("pgumbel", q, location, scale, lower.tail, max,
          PACKAGE="Rufus")

##' @param log logical; if \code{TRUE}, probabilities p are given as
##' log(p).
##' @rdname pgumbel
##' @export
dgumbel <-
  function(x, location = 0, scale = 1, log = FALSE, max = TRUE)
    .Call("dgumbel", x, location, scale, log, max, PACKAGE="Rufus")

##' @param n number of observations.
##' @rdname pgumbel
##' @export
rgumbel <-
    function(n, location = 0, scale = 1, max = TRUE) {
  if(max)
    location - scale * log(-log(runif(n)))
  else
    location + scale * log(-log(runif(n)))
}

##' @param p vector of probabilities.
##' @rdname pgumbel
##' @export
qgumbel <-
  function(p, location = 0, scale = 1, lower.tail = TRUE, max = TRUE)
{
  if(max)  ## right skew, loglog link
    location - scale * log(-log(p))
  else ## left skew, cloglog link
    location + scale * log(-log(1 - p))
}

##' @rdname pgumbel
##' @export
ggumbel <- function(x, max = TRUE)
    .Call("ggumbel", x, max, PACKAGE="Rufus")
