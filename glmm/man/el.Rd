\name{el}
\alias{el}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
el(Y, X, eta, family.mcml)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
%%     ~~Describe \code{Y} here~~
}
  \item{X}{
%%     ~~Describe \code{X} here~~
}
  \item{eta}{
%%     ~~Describe \code{eta} here~~
}
  \item{family.mcml}{
%%     ~~Describe \code{family.mcml} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

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
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
