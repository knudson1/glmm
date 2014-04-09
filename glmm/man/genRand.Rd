\name{genRand}
\alias{genRand}

\title{Random effect generator for Monte Carlo likelihood approximation}
\description{Generates m vectors of random effects from a normal distribution based on PQL estimates for the mean and variance.}
\usage{genRand(sigma.star, s.star, z.list, m)}

\arguments{
  \item{sigma.star}{A vector of distinct standard deviations. Length of sigma.star is the same as the length of z.list.}
  \item{s.star}{A vector containing the means.}
  \item{z.list}{A list of model matrices for the random effects. Length of z.list is the same as the length of sigma.star.}
  \item{m}{The Monte Carlo sample size, which is the number of vectors to be drawn from the specified normal distribution.}
}
\details{Generates a vector of random effects from a standard normal distribution,  centers the vector at the PQL mean u.star, and scales the vector using the PQL variance components.
}
\value{
\item{u}{A matrix of random effects. Each row of the matrix is one vector of random effects generated from the normal distribution specified.}
\item{u.star}{The mean vector of the generated random effects.}
 }
\references{
%% ~put references to the literature/web site here ~
}
\author{Christina Knudson}
\note{ Not intended for use by naive users.  Use \code{\link{mcml}},
  which calls them.}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{mcml}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (sigma.star, s.star, z.list, m) 
{
    nrand <- lapply(z.list, ncol)
    nrandom <- unlist(nrand)
    q <- sum(nrandom)
    if (q != length(s.star)) 
        stop("Number of random effects in dispute between s.star and mod.mcml$z")
    eek <- getEk(z.list)
    Aks <- Map("*", eek, sigma.star)
    A.star <- addVecs(Aks)
    u.star <- A.star * s.star
    u <- matrix(data = NA, nrow = m, ncol = q)
    for (k in 1:m) {
        u.swoop <- rnorm(q)
        u[k, ] <- u.swoop * A.star + u.star
    }
    list(u = u, u.star = u.star)
  }
}
\keyword{ generalized linear mixed model }
\keyword{Monte Carlo likelihood approximation }
