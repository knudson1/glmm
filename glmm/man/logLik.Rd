\name{logLik.glmm}
\alias{logLik.glmm}
\concept{generalized linear mixed model}
\concept{Monte Carlo}
\concept{likelihood approximation}

\title{Monte Carlo Log Likelihood}
\description{A function that calculates the Monte Carlo log likelihood evaluated at the the Monte Carlo maximum likelihood estimates returned from \code{\link{glmm}}.

}
\usage{
\method{logLik}{glmm}(object,\dots)
}

\arguments{
  \item{object}{An object of class \code{glmm} usually created using \code{\link{glmm}}.}
  \item{\dots}{further arguments passed to or from other methods.}
}

\value{
 \item{logLik }{The variance-covariance matrix for the parameter estimates}
}

\author{Christina Knudson}

\seealso{
\code{\link{glmm}} for model fitting.
}
\examples{
library(glmm)
data(BoothHobert)
set.seed(1234)
mod <- glmm(y~0+x1, list(y~0+z1), varcomps.names=c("z1"), 
data=BoothHobert, family.glmm=bernoulli.glmm, m=100)
logLik(mod)
}

\keyword{models}
