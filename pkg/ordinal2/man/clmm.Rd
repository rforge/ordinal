\name{clmm}
\alias{clmm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Cumulative Link Mixed Models
}
\description{
  Fits Cumulative Link Mixed Models with one or more random effects via
  the Laplace approximation
}
\usage{
clmm(formula, data, weights, start, subset, na.action, contrasts, Hess =
TRUE, model = TRUE, link = c("logit", "probit", "cloglog", "loglog",
"cauchit"), doFit = TRUE, control, threshold = c("flexible",
"symmetric", "equidistant"), ...)

%% also document getNLA(rho, par) here and include examples
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
    a two-sided linear formula object describing the fixed-effects part
    of the model, with the response on the left of a ~ operator and the
    terms, separated by + operators, on the right. The vertical bar
    character "|" separates an expression for a model matrix and a
    grouping factor. 
  }
  \item{data}{
    an optional data frame in which to interpret the variables occurring
    in the formula.
  }
  \item{weights}{
    optional case weights in fitting. Defaults to 1.
  }
  \item{start}{
    optional initial values for the parameters in the format
    \code{c(alpha, beta, tau)}, where \code{alpha} are the threshold
    parameters, \code{beta} are the fixed regression parameters and
    \code{tau} are variance parameters for the random effects on the log
    scale.
  }
  \item{subset}{
    expression saying which subset of the rows of the data should  be
    used in the fit. All observations are included by default. 
  }
  \item{na.action}{
    a function to filter missing data.
  }
  \item{contrasts}{
    a list of contrasts to be used for some or all of
    the factors appearing as variables in the model formula.
  }
  \item{Hess}{
    logical for whether the Hessian (the inverse of the observed
    information matrix)
    should be computed.
    Use \code{Hess = TRUE} if you intend to call \code{summary} or
    \code{vcov} on the fit and \code{Hess = FALSE} in all other instances
    to save computing time.
  }
  \item{model}{
    logical for whether the model frames should be part of the returned
    object.
  }
  \item{link}{
    link function, i.e. the type of location-scale distribution
    assumed for the latent distribution. The default \code{"logit"} link
    gives the proportional odds mixed model.
  }
  \item{doFit}{
    logical for whether the model should be fit or the model
    environment should be returned.
  }
  \item{control}{
    a call to \code{\link{clmm.control}}
  }
  \item{threshold}{
    specifies a potential structure for the thresholds
    (cut-points). \code{"flexible"} provides the standard unstructured
    thresholds, \code{"symmetric"} restricts the distance between the
    thresholds to be symmetric around the central one or two thresholds
    for odd or equal numbers or thresholds respectively, and
    \code{"equidistant"} restricts the distance between consecutive
    thresholds to the same value.
  }
  \item{\dots}{
    additional arguments are passed on to \code{\link{clm.control}}.
  }
}
\details{
  There are standard print and summary methods implemented for for
  \code{"clmm"} objects.
}
\value{ a list containing
  \item{alpha}{threshold parameters.}
  \item{beta}{fixed effect regression parameters.}
  \item{stDev}{standard deviation of the random effect terms.}
  \item{tau}{\code{log(stDev)} - the scale at which the log-likelihood
    function is optimized.}
  \item{coefficients}{the estimated model parameters = \code{c(alpha,
      beta, tau}.}
  \item{Hessian}{Hessian of the model coefficients.}
  \item{edf}{the estimated degrees of freedom used by the model =
    \code{length(coefficients)}.} 
  \item{nobs}{\code{sum(weights)}.}
  \item{n}{length(y).}
  \item{fitted.values}{fitted values evaluated with the random effects
    at their conditional modes.}
  \item{df.residual}{residual degrees of freedom; \code{length(y) -
      sum(weights)}} 
  \item{tJac}{Jacobian of the threshold function corresponding to the
    mapping from standard flexible thresholds to those used in the
    model.}
  \item{terms}{the terms object for the fixed effects.}
  \item{contrasts}{contrasts applied to the fixed model terms.}
  \item{na.action}{the function used to filter missing data.}
  \item{call}{the matched call.}
  \item{logLik}{value of the log-likelihood function for the model at
    the optimum.}
  \item{Niter}{number of Newton iterations in the inner loop update of
    the conditional modes of the random effects.}
  \item{optRes}{list of results from the optimizer.}
  \item{ranef}{list of the conditional modes of the random effects.}
  \item{condVar}{list of the conditional variance of the random effects
    at their conditional modes.}
  \item{}{}
  \item{}{}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{ bla
%% ~put references to the literature/web site here ~
}
\author{
  Rune Haubo B Christensen
}
\examples{

## Get test data:
data(soup, package = "ordinal")

## Cumulative link mixed model with two random terms:
mm1 <- clmm(SURENESS ~ PROD + (1|RESP) + (1|RESP:PROD), data = soup,
            link = "probit", threshold = "equidistant")
mm1
summary(mm1)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}

