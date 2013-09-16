\name{gnorm}
\alias{gcauchy}
\alias{gfun}
\alias{glogis}
\alias{gnorm}
\title{Gradients of common densities}
\usage{
  gnorm(x, mean = 0, sd = 1)

  glogis(x)

  gcauchy(x)
}
\arguments{
  \item{x}{numeric vector of quantiles.}

  \item{mean}{numeric scalar; mean.}

  \item{sd}{numeric scalar; standard deviation.}
}
\value{
  a numeric vector of gradients.
}
\description{
  Gradients of common density functions in their standard
  forms, i.e., with zero location (mean) and unit scale.
  These are implemented in C++ for speed and care is taken
  that the correct results are provided for the argument
  being \code{NA}, \code{NaN}, \code{Inf}, \code{-Inf} or
  just extremely small or large.
}
\details{
  The gradients are given by: \itemize{ \item{gnorm: If
  \eqn{f(x)} is the normal density with mean 0 and spread
  1, then the gradient is \deqn{f'(x) = -x f(x)} }
  \item{glogis: If \eqn{f(x)} is the logistic density with
  mean 0 and scale 1, then the gradient is \deqn{f'(x) = 2
  \exp(-x)^2 (1 + \exp(-x))^{-3} -
  \exp(-x)(1+\exp(-x))^{-2}}{f'(x) = 2 exp(-x)^2 (1 +
  exp(-x))^{-3} - exp(-x)(1+exp(-x))^{-2}} } \item{pcauchy:
  If \eqn{f(x) = [\pi(1 + x^2)^2]^{-1}}{f(x) =1 / [pi (1 +
  x^2)^2]} is the cauchy density with mean 0 and scale 1,
  then the gradient is \deqn{f'(x) = -2x [\pi(1 +
  x^2)^2]^{-1}}{f'(x) = -2x / [pi (1 + x^2)^2]} } }
}
\examples{
x <- -5:5
gnorm(x)
glogis(x)
gcauchy(x)

## Change default mean and standard deviation:
gnorm(x, 1, 2)
}
\author{
  Rune Haubo Bojesen Christensen
}
\seealso{
  Gradients of densities are also implemented for the
  extreme value distribtion
  (\code{\link[=dgumbel]{gumbel}}) and the the log-gamma
  distribution (\code{\link[=lgamma]{log-gamma}}).

  Other distributions: \code{\link{dgumbel}},
  \code{\link{ggumbel}}, \code{\link{pgumbel}},
  \code{\link{qgumbel}}, \code{\link{rgumbel}}
}
\keyword{distribution}

