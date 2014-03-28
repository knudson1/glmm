\name{getFamily}
\alias{getFamily}

\title{Determines the MCML family.}
\description{Determines the MCML family specified by the user.}
\usage{getFamily(family.mcml)}

\arguments{
  \item{family.mcml}{The MCML family. Must be of class mcml.family.}
}
\details{The input family.mcml can either by a string, a function, or the list of functions contained in family.mcl. This function determines the MCML family and calls those functions.}
\value{Returns the list of functions in the specified MCML family.}
\references{
%% ~put references to the literature/web site here ~
}
\author{Christina Knudson}
\note{ Not intended for use by naive users.  Use \code{\link{mcml}},
  which calls them.}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{mcml}}
 \code{\link{bernoulli.mcml}}
 \code{\link{poisson.mcml}}
}
\examples{
getFamily(bernoulli.mcml)
getFamily(poisson.mcml)
## The function is currently defined as
function (family.mcml) 
{
    if (is.character(family.mcml)) 
        family.mcml <- get(family.mcml, mode = "function", envir = parent.frame())
    if (is.function(family.mcml)) 
        family.mcml <- family.mcml()
    if (!inherits(family.mcml, "mcml.family")) 
        stop(" 'family' not recognized")
    family.mcml
  }
}

\keyword{ generalized linear mixed model }
\keyword{Monte Carlo likelihood approximation }
