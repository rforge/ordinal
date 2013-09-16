\name{pgumbel}
\alias{dgumbel}
\alias{ggumbel}
\alias{ggumbel,}
\alias{gumbel}
\alias{pgumbel}
\alias{qgumbel}
\alias{rgumbel}
\title{The Gumbel distribution}
\usage{
  pgumbel(q, location = 0, scale = 1, lower.tail = TRUE,
    max = TRUE)

  dgumbel(x, location = 0, scale = 1, log = FALSE,
    max = TRUE)

  rgumbel(n, location = 0, scale = 1, max = TRUE)

  qgumbel(p, location = 0, scale = 1, lower.tail = TRUE,
    max = TRUE)

  ggumbel(x, max = TRUE)
}
\arguments{
  \item{q,x}{numeric vector of quantiles.}

  \item{location}{numeric scalar.}

  \item{scale}{numeric scalar.}

  \item{lower.tail}{logical; if \code{TRUE} (default),
  probabilities are \eqn{P[X \leq x]}{P[X <= x]} otherwise,
  \eqn{P[X > x]}.}

  \item{max}{distribution for extreme maxima (default) or
  minima? The default corresponds to the standard
  right-skew Gumbel distribution.}

  \item{log}{logical; if \code{TRUE}, probabilities p are
  given as log(p).}

  \item{n}{number of observations.}

  \item{p}{vector of probabilities.}
}
\value{
  \code{pgumbel} gives the distribution function,
  \code{dgumbel} gives the density, \code{ggumbel} gives
  the gradient of the density, \code{qgumbel} is the
  quantile function, and \code{rgumbel} generates random
  deviates.
}
\description{
  Density, distribution function, quantile function, random
  generation, and gradient of density of the extreme value
  (maximum and minimum) distributions. The Gumbel
  distribution is also known as the extreme value maximum
  distribution, the double-exponential distribution and the
  log-Weibull distribution.
}
\details{
  \code{dgumbel}, \code{pgumbel} and \code{ggumbel} are
  implemented in C++ for speed and care is taken that
  'correct' results are provided for values of \code{NA},
  \code{NaN}, \code{Inf}, \code{-Inf} or just extremely
  small or large.
}
\examples{
## Illustrating the symmetry of the distribution functions:
pgumbel(5) == 1 - pgumbel(-5, max=FALSE) ## TRUE
dgumbel(5) == dgumbel(-5, max=FALSE) ## TRUE
ggumbel(5) == -ggumbel(-5, max=FALSE) ## TRUE

## More examples:
x <- -5:5

(pp <- pgumbel(x))
qgumbel(pp)
dgumbel(x)
ggumbel(x)

(ppp <- pgumbel(x, max=FALSE))
## Observe lack of precision:
qgumbel(ppp, max=FALSE)
dgumbel(x, max=FALSE)
ggumbel(x, max=FALSE)

## random deviates:
set.seed(1)
(r1 <- rgumbel(10))
set.seed(1)
r2 <- -rgumbel(10, max = FALSE)
all(r1 == r2) ## TRUE
}
\author{
  Rune Haubo Bojesen Christensen
}
\references{
  \url{wikipedia.org/wiki/Gumbel_distribution}
}
\seealso{
  Gradients of densities are also implemented for the
  normal, logistic, cauchy, cf. \code{\link[=gnorm]{gfun}}
  and the log-gamma distribution, cf. \code{\link{lgamma}}.

  Other distributions: \code{\link{gcauchy}},
  \code{\link{glogis}}, \code{\link{gnorm}}
}
\keyword{distribution}

