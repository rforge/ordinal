\name{clm.fit}
\alias{clm.fit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  fit cumulative link functions
%%  ~~function to do ... ~~
}
\description{
  A direct fitter
}
\usage{
clm.fit(y, X, weights = rep(1, nrow(X)), offset = rep(0, nrow(X)),
control = list(), start,
link = c("logit", "probit", "cloglog", "loglog", "cauchit"),
threshold = c("flexible", "symmetric", "equidistant"))

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{X}{
%%     ~~Describe \code{X} here~~
}
  \item{weights}{
%%     ~~Describe \code{weights} here~~
}
  \item{offset}{
%%     ~~Describe \code{offset} here~~
}
  \item{control}{
%%     ~~Describe \code{control} here~~
}
  \item{start}{
%%     ~~Describe \code{start} here~~
}
  \item{link}{
%%     ~~Describe \code{link} here~~
}
  \item{threshold}{
%%     ~~Describe \code{threshold} here~~
}
}
\details{ bla
%%  ~~ If necessary, more details than the description above ~~
}
\value{ bla
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
%%  ~~who you are~~
}
\note{ bla
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ bla
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
