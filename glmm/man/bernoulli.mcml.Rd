\name{bernoulli.mcml}
\alias{bernoulli.mcml}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Functions for the bernoulli family. }
\description{
Given a scalar  eta, this calculates the cumulant and two derivatives for the bernoulli family. Also checks that the data are entered correctly.}
\usage{
bernoulli.mcml()
}
%- maybe also 'usage' for other objects documented here.
\details{
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
function () 
{
    family.mcml <- "bernoulli.mcml"
    cum <- function(eta) log(1 + exp(eta))
    cp <- function(eta) exp(eta)/(1 + exp(eta))
    cpp <- function(eta) exp(eta)/(1 + exp(eta)) - exp(2 * eta) * 
        (1 + exp(eta))^(-2)
    checkData <- function(x) {
        bads <- sum(x > 1) + sum(x < 0)
        if (bads > 0) 
            stop("data must be between 0 and 1.")
        return(NULL)
    }
    out <- list(family.mcml = family.mcml, cum = cum, cp = cp, 
        cpp = cpp, checkData = checkData)
    class(out) <- "mcml.family"
    return(out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ generalized linear mixed model }
\keyword{ Monte Carlo likelihood approximation }% 