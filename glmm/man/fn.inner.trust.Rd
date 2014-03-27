\name{fn.inner.trust}
\alias{fn.inner.trust}
\title{The penalized likelihood}
\description{Penalized quasi-likelihood (PQL) consists of an inner optimization and an outer optimization. This is the inner function, which is a penalized likelihood.  This is to be optimized by trust.}
\usage{
fn.inner.trust(mypar, Y, X, Z, A, family.mcml, nbeta)
}
\arguments{
  \item{mypar}{A vector of parameters to be optimized: the fixed effects (beta) and the standardized random effects (s).}
  \item{Y}{The response vector.}
  \item{X}{The model matrix for the fixed effects.}
  \item{Z}{The model matrix for the random effects.}
  \item{A}{A matrix containing the standard deviations (the square root of the variance components nu).}
  \item{family.mcml}{The name of the family. Must be class mcml.family}
  \item{beta}{The length of the beta vector.}
}
\details{Evaluates the penalized log likelihood and two derivatives.}
\value{
\item{value }{The value of the penalized log likelihood}
\item{gradient }{The gradient vector of the  penalized log likelihood}
\item{hessian }{The hessian matrix of the  penalized log likelihood}
}
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
function (mypar, Y, X, Z, A, family.mcml, nbeta) 
{
    beta <- mypar[1:nbeta]
    s <- mypar[-(1:nbeta)]
    eta <- X \%*\% beta + Z \%*\% A \%*\% s
    family.mcml <- getFamily(family.mcml)
    mu <- family.mcml$cp(eta)
    value <- el(Y, X, eta, family.mcml)$value - 0.5 * s \%*\% s
    db <- t(X) \%*\% (Y - mu)
    ds <- A \%*\% t(Z) \%*\% (Y - mu) - s
    gradient <- c(db, ds)
    cdub <- family.mcml$cpp(eta)
    cdub <- as.vector(cdub)
    cdub <- diag(cdub)
    kyoo <- nrow(A)
    piece1 <- (-t(X) \%*\% cdub \%*\% X)
    piece2 <- (-t(X) \%*\% cdub \%*\% Z \%*\% A)
    piece3 <- (-A \%*\% t(Z) \%*\% cdub \%*\% X)
    piece4 <- (-A \%*\% t(Z) \%*\% cdub \%*\% Z \%*\% A - diag(kyoo))
    hessian <- rbind(cbind(piece1, piece2), cbind(piece3, piece4))
    list(value = value, gradient = gradient, hessian = hessian)
  }
}

\keyword{ generalized linear mixed model }
