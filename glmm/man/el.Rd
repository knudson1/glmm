\name{el}
\alias{el}
\title{Log likelihood of the data}
\description{Calculates the log density of the data, its gradient vector, and its Hessian matrix.}
\usage{el(Y, X, eta, family.mcml)}
\arguments{
  \item{Y}{Response vector}
  \item{X}{Model matrix for the fixed effects}
  \item{eta}{Canonical parameter vector}
  \item{family.mcml}{The name of the family. Must be class mcml.family}
}
\details{The lengths of Y and eta should be equal.}
\value{
\item{value }{The value of the log density of the data.}
\item{gradient }{The gradient vector of the log density of the data.}
\item{hessian }{The hessian matrix of the log density of the data.}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{Christina Knudson}
\note{}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (Y, X, eta, family.mcml) 
{
    family.mcml <- getFamily(family.mcml)
    value <- Y \%*\% eta - sum(family.mcml$cum(eta))
    mu <- family.mcml$cp(eta)
    gradient <- t(X) \%*\% (Y - mu)
    cdub <- as.vector(family.mcml$cpp(eta))
    cdubmat <- diag(cdub)
    hessian <- t(X) \%*\% (-cdubmat) \%*\% X
    list(value = value, gradient = gradient, hessian = hessian)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ generalized linear mixed model}
\keyword{ monte carlo likelihood approximation }
