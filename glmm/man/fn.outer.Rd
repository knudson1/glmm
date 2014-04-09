\name{fn.outer}
\alias{fn.outer}
\title{The penalized quasi-likelihood}
\description{Penalized quasi-likelihood (PQL) consists of an inner optimization and an outer optimization. This is the outer function, the penalized quasi-likelihood. This is to be optimized using Nelder-Mead.}
\usage{fn.outer(par, beta, s, Y, X, Z, eek, family.mcml)}
\arguments{
  \item{par}{The vector of standard deviations (square root of the variance components nu).}
  \item{beta}{The vector of fixed effects.}
  \item{s}{The vector of random effects.}
  \item{Y}{The response vector.}
  \item{X}{The model matrix for the fixed effects.}
  \item{Z}{The model matrix for the random effects.}
  \item{eek}{A list of vectors containing 0 or 1.}
  \item{family.mcml}{The name of the family. Must be class mcml.family}
}
\details{}
\value{The evaluated penalized quasi-likelihood.}
\references{
%% ~put references to the literature/web site here ~
}
\author{Christina Knudson}
\note{ Not intended for use by naive users.  Use \code{\link{mcml}},
  which calls them.}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{fn.outer}}
 \code{\link{mcml}}
 \code{\link{pql}}

}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (par, beta, s, Y, X, Z, eek, family.mcml) 
{
    sigma <- par
    Aks <- Map("*", eek, sigma)
    A <- addVecs(Aks)
    A <- diag(A)
    nbeta <- length(beta)
    inner.optim <- trust(fn.inner.trust, parinit = c(beta, s), 
        rinit = 5, rmax = 10000, minimize = F, Y = Y, X = X, 
        Z = Z, A = A, nbeta = nbeta, family.mcml = family.mcml)
    parms <- inner.optim$argument
    beta.twid <<- parms[1:nbeta]
    s.twid <<- parms[-(1:nbeta)]
    eta.twid <- (X \%*\% beta.twid + Z \%*\% A \%*\% s.twid)
    family.mcml <- getFamily(family.mcml)
    piece1 <- el(Y, X, eta.twid, family.mcml)$value
    dubya <- family.mcml$cpp(eta.twid)
    W <- diag(as.vector(dubya))
    thatthing <- A \%*\% t(Z) \%*\% W \%*\% Z \%*\% A + diag(nrow(A))
    L <- chol(thatthing)
    piece2a <- sum(log(diag(L)))
    piece2b <- 0.5 * s \%*\% s
    minvalue <- piece2a + piece2b - piece1
    as.vector(minvalue)
  }
}
\keyword{ generalized linear mixed model }

