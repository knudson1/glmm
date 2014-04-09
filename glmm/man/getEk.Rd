\name{getEk}
\alias{getEk}
\title{Function relating random effects and variance components.}
\description{Determines the relationships between random effects and variance components. }
\usage{getEk(z)}

\arguments{
 \item{z}{A list of design matrices for the random effects.}
}
\details{Determines the relationships between random effects and variance components. If the ith random effect has variance component j, then the jth vector will have a 1 in the ith position. Thus the vectors' elements are either 0 or 1. The sum of the vectors is a vector containing only 1s.}
\value{A list of vectors with components either 0 or 1.}
\references{
%% ~put references to the literature/web site here ~
}
\author{Christina Knudson}
\note{ Not intended for use by naive users.  Use \code{\link{mcml}},
  which calls them.}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{mcml}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (z) 
{
    nrand <- lapply(z, ncol)
    nrandom <- unlist(nrand)
    totnrandom <- sum(nrandom)
    nsigma <- length(z)
    Ekmat <- matrix(data = c(0), nrow = nsigma, ncol = totnrandom)
    Ek <- list()
    cumnrand <- c()
    cumnrand[1] <- 0
    for (i in 1:length(nrandom)) {
        cumnrand[i + 1] <- sum(nrandom[1:i])
    }
    for (i in 1:nsigma) {
        for (j in 1:totnrandom) {
            if (cumnrand[i] < j & cumnrand[i + 1] > j - 1) 
                Ekmat[i, j] <- 1
        }
        Ek[[i]] <- Ekmat[i, ]
    }
    Ek
  }
}
\keyword{ generalized linear mixed model }
\keyword{Monte Carlo likelihood approximation }
