\name{bernoulli.mcml}
\alias{bernoulli.mcml}

\title{
Functions for the Bernoulli family. }
\description{
Given a scalar eta, this calculates the cumulant and two derivatives for the Bernoulli family. Also checks that the data are entered correctly.}
\usage{
bernoulli.mcml()
}

\details{
 }
\value{
\item{family.mcml }{The family name, as a string.}
  \item{cum }{The cumulant function.}
  \item{cp }{The first derivative of the cumulant function.}
  \item{cpp }{The second derivative of the cumulant function.}
  \item{checkData }{A function to check that all data are either 0 or 1.}
}
\references{}
\author{Christina Knudson}
\note{This function is to be used by the mcml command.}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{mcml}}
}
\examples{
bernoulli.mcml()$family.mcml
bernoulli.mcml()$cum(2)
bernoulli.mcml()$cp(2)
bernoulli.mcml()$cpp(2)

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

\keyword{ generalized linear mixed model }
\keyword{ Monte Carlo likelihood approximation }
