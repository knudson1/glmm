\name{distRand}
\alias{distRand}
\title{Distribution of the random effects}
\description{Random effects are assumed normally distributed. This function calculates the log density of the random effects, its gradient vector, and its Hessian matrix.}
\usage{distRand(nu, U, z.list, mu)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nu}{A vector of the distinct variance components.}
  \item{U}{A vector of (generated) random effects.}
  \item{z.list}{A list of the model matrices for the random effects.}
  \item{mu}{The mean vector for the random effects.}
}
\details{The arguments mu and U should be vectors of equal length. The length of z.list should be equal to the length of nu.}

\value{
\item{value }{The value of the log density of the random effects.}
\item{gradient }{The gradient vector of the log density of the random effects.}
\item{hessian }{The hessian matrix of the log density of the random effects.}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{Christina Knudson}
\note{  Not intended for use by naive users.  Use \code{\link{mcml}},
  which calls them.}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (nu, U, z.list, mu) 
{
    T <- length(z.list)
    nrand <- lapply(z.list, ncol)
    nrandom <- unlist(nrand)
    totnrandom <- sum(nrandom)
    mu.list <- U.list <- NULL
    if (T == 1) {
        U.list[[1]] <- U
        mu.list[[1]] <- mu
    }
    if (T > 1) {
        U.list[[1]] <- U[1:nrandom[1]]
        mu.list[[1]] <- mu[1:nrandom[1]]
        for (t in 2:T) {
            thing1 <- sum(nrandom[1:t - 1]) + 1
            thing2 <- sum(nrandom[1:t])
            U.list[[t]] <- U[thing1:thing2]
            mu.list[[t]] <- mu[thing1:thing2]
        }
    }
    val <- gradient <- Hessian <- rep(0, T)
    for (t in 1:T) {
        you <- as.vector(U.list[[t]])
        mew <- as.vector(mu.list[[t]])
        Umu <- (you - mew) \%*\% (you - mew)
        val[t] <- as.numeric(-0.5 * nrandom[t] * log(nu[t]) - 
            Umu/(2 * nu[t]))
        gradient[t] <- -nrandom[t]/(2 * nu[t]) + Umu/(2 * (nu[t])^2)
        Hessian[t] <- nrandom[t]/(2 * (nu[t])^2) - Umu/((nu[t])^3)
    }
    value <- sum(val)
    if (T > 1) 
        hessian <- diag(Hessian)
    if (T == 1) 
        hessian <- matrix(Hessian, nrow = 1, ncol = 1)
    list(value = value, gradient = gradient, hessian = hessian)
  }
}

\keyword{ generalized linear mixed model}
\keyword{ monte carlo likelihood approximation }
