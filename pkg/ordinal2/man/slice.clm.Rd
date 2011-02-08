\name{slice.clm}
\alias{slice.clm}
\alias{plot.slice.clm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  slice the likelihood of a clm
}
\description{
  slice likelihood and plot the slice
}
\usage{
slice.clm(x, which = seq_along(par), lambda = 3, grid = 100,
quad = TRUE)

plot.slice.clm(x, which = seq_along(x), type = c("quadratic", "linear"),
plot.mle = TRUE, ...) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    %%     ~~Describe \code{x} here~~
  }
  \item{which}{
    %%     ~~Describe \code{which} here~~
  }
  \item{lambda}{
    %%     ~~Describe \code{lambda} here~~
  }
  \item{grid}{
    %%     ~~Describe \code{grid} here~~
  }
  \item{quad}{
    %%     ~~Describe \code{quad} here~~
  }
  \item{type}{
    %%     ~~Describe \code{type} here~~
  }
  \item{plot.mle}{
    %%     ~~Describe \code{plot.mle} here~~
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
\author{
  Rune Haubo B Christensen
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}

