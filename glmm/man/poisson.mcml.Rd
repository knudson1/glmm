\name{poisson.mcml}
\alias{poisson.mcml}

\title{
Functions for the Poisson family.
}
\description{
Given a scalar eta, this calculates the cumulant and two derivatives for the Poisson family. Also checks that the data are entered correctly.
}
\usage{
poisson.mcml()
}
%- maybe also 'usage' for other objects documented here.
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{family.mcml }{The family name, as a string.}
  \item{cum }{The cumulant function.}
  \item{cp }{The first derivative of the cumulant function.}
  \item{cpp }{The second derivative of the cumulant function.}
  \item{checkData }{A function to check that all data are either 0 or 1.}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Christina Knudson
}
\note{
This function is typically called through the getFamily function.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{getFamily}}
}
\examples{
poisson.mcml()$family.mcml
poisson.mcml()$cum(2)
poisson.mcml()$cp(2)
poisson.mcml()$cpp(2)

## The function is currently defined as
function () 
{
    family.mcml <- "poisson.mcml"
    cum <- function(eta) exp(eta)
    cp <- function(eta) exp(eta)
    cpp <- function(eta) exp(eta)
    checkData <- function(x) {
        bads <- sum(x < 0)
        if (bads > 0) 
            stop("data must be nonnegative.")
        return(NULL)
    }
    out <- list(family.mcml = family.mcml, cum = cum, cp = cp, 
        cpp = cpp, checkData = checkData)
    class(out) <- "mcml.family"
    return(out)
  }
}

\keyword{ generalized linear mixed model }
\keyword{ Monte Carlo likelihood approximation }
