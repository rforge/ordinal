\name{clm}
\alias{clm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Cumulative Link Functions
%%  ~~function to do ... ~~
}
\description{
  Fits cumulative link models (CLMs) such as the propotional odds
  model. The model allows for various link functions and structured
  thresholds that restricts the thresholds or cut-points to be e.g.,
  equidistant or symmetrically arranged around the central
  threshold(s). A Newton algorithm is used to optimize the likelihood
  function which is very speedy compared to general optimization
  routines such as \code{\link{optim}}.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
clm(formula, data, weights, start, subset, doFit = TRUE, na.action,
  contrasts, model = TRUE, control, link = c("logit", "probit",
  "cloglog", "loglog", "cauchit"), threshold = c("flexible",
  "symmetric", "equidistant"), ...) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
    a formula expression as for regression models, of the form
    \code{response ~ predictors}. The response should be a factor
    (preferably an ordered factor), which will be interpreted as an
    ordinal response with levels ordered as in the factor.
    The model must have an intercept: attempts to remove one will
    lead to a warning and will be ignored. An offset may be used.  See the
    documentation of \code{\link{formula}} for other details.
  }
  \item{data}{
    an optional data frame in which to interpret the variables occurring
    in the formulas.
  }
  \item{weights}{
    optional case weights in fitting. Defaults to 1.
  }
  \item{start}{
    initial values for the parameters in the format
    \code{c(alpha, beta)}, where \code{alpha} are the threshold
    parameters and \code{beta} are the regression parameters.
  }
  \item{subset}{
    expression saying which subset of the rows of the data should  be used
    in the fit. All observations are included by default.
  }
  \item{doFit}{
    logical for whether the model should be fit or the model
    environment should be returned.
  }
  \item{na.action}{
    a function to filter missing data. Applies to terms in all three formulae.
  }
  \item{contrasts}{
    a list of contrasts to be used for some or all of
    the factors appearing as variables in the model formula.
  }
  \item{model}{
    logical for whether the model frames should be part of the returned
    object.
  }
  \item{control}{
    a call to \code{\link{clm.control}}.
  }
  \item{link}{
    link function, i.e. the type of location-scale distribution
    assumed for the latent distribution. The default \code{"logit"} link
    gives the proportional odds model.
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
  There are methods for the standard model-fitting functions, including
  \code{\link{summary}},
  \code{\link{logLik}},
  \code{\link{vcov}} and
  \code{\link[=slice.clm]{slice}}.
  %% \code{slice}.

  
}
\value{
  A list containing ...
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\author{
  Rune Haubo B Christensen
}
\keyword{models}

