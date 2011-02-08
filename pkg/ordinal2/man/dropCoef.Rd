\name{dropCoef}
\alias{dropCoef}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Insure Full Rank Design Matrix
}
\description{
  Coefficients (columns) are removed from a design matrix, equivalently
  \code{model.matrix} to insure that it has full rank.
}
\usage{
dropCoef(X)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
    a design matrix, e.g., the result of \code{\link{model.matrix}}
    possibly of less than full rank, i.e., with redundant parameters
    (columns). 
  }
}
\details{
  The redundant columns of the design matrix are identified with the
  LAPACK implementation of the \code{\link{qr}} decomposition. 
}
\value{
  the design matrix \code{X} without redundant columns.
}
\author{
  Rune Haubo B Christensen
}
\seealso{
  \code{\link{qr}} and \code{\link{lm}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}

