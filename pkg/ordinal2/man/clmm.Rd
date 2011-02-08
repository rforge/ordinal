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
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
%%     ~~Describe \code{formula} here~~
}
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{weights}{
%%     ~~Describe \code{weights} here~~
}
  \item{start}{
%%     ~~Describe \code{start} here~~
}
  \item{subset}{
%%     ~~Describe \code{subset} here~~
}
  \item{na.action}{
%%     ~~Describe \code{na.action} here~~
}
  \item{contrasts}{
%%     ~~Describe \code{contrasts} here~~
}
  \item{Hess}{
%%     ~~Describe \code{Hess} here~~
}
  \item{model}{
%%     ~~Describe \code{model} here~~
}
  \item{link}{
%%     ~~Describe \code{link} here~~
}
  \item{doFit}{
%%     ~~Describe \code{doFit} here~~
}
  \item{control}{
%%     ~~Describe \code{control} here~~
}
  \item{threshold}{
%%     ~~Describe \code{threshold} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
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

