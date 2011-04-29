\name{gumbel}
\alias{dgumbel}
\alias{pgumbel}
\alias{ggumbel}
\alias{dgumbel2}
\alias{pgumbel2}
\alias{ggumbel2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  The Gumbel distribution
%%  ~~function to do ... ~~
}
\description{
  Density, distribution function and gradient of density of the extreme
  value (minimum and maximum) distributions. The extreme value maximum
  distribution is also known as the Gumbel distribution.
  These are implemented in C
  for speed and care is taken that the correct results are provided for
  values of \code{NA}, \code{NaN}, \code{Inf}, \code{-Inf} or just
  extremely small or large. 
  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{

pgumbel(q, location = 0, scale = 1, lower.tail = TRUE)

pgumbel2(q, location = 0, scale = 1, lower.tail = TRUE)

dgumbel(x, location = 0, scale = 1, log = FALSE)

dgumbel2(x, location = 0, scale = 1, log = FALSE)

ggumbel(x)

ggumbel2(x)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x,q}{
    numeric vector of quantiles.
  }
  \item{location}{
    numeric scalar.
  }
  \item{scale}{
    numeric scalar.
  }
  \item{lower.tail}{
    logical; if \code{TRUE} (default), probabilities are
    \eqn{P[X \leq x]}{P[X <= x]} otherwise, \eqn{P[X > x]}.
  }
  \item{log}{
    logical; if \code{TRUE}, probabilities p are given as log(p).
  }
}
\details{

  The Gumbel (maximum) distribution function is defined as \ldots and
  corresponds to the complementary(?) log-log (\code{cloglog}) link.

  The density and gradient of the density are defined as\ldots
  
  These distribution functions, densities and gradients are used in the
  Newton-Raphson algorithms in fitting cumulative link models with
  \code{\link{clm}} and cumulative link mixed models with
  \code{\link{clmm}}. 
}
\value{
  \code{pgumbel}[2] gives the distribution function, \code{dgumbel}[2]
  gives the density and \code{ggumbel}[2] gives the gradient of the
  density. 
}
\references{
  \url{wikipedia.org/wiki/Gumbel_distribution}
}
\seealso{
  Gradients of densities are also implemented for the normal, logistic,
  cauchy, log-gamma distribution and Aranda-Ordaz
  distributions. 
}
\author{
  Rune Haubo B Christensen
}
\examples{

## Illustrating the symmetry of the distribution functions:
pgumbel(5) == 1 - pgumbel2(-5) ## TRUE
dgumbel(5) == dgumbel2(-5)     ## TRUE
ggumbel(5) == -ggumbel2(-5)    ## TRUE

## More examples:
x <- -5:5

pgumbel(x)
dgumbel(x)
ggumbel(x)

pgumbel2(x)
dgumbel2(x)
ggumbel2(x)

}
\keyword{models}

