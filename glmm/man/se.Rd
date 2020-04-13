\name{se}
\alias{se}
\concept{generalized linear mixed model}

\title{Standard Error}
\description{A function that calculates the standard error for the Monte Carlo maximum likelihood estimates returned from \code{\link{glmm}}.

}
\usage{se(object)}

\arguments{
  \item{object}{An object of class \code{glmm} usually created using \code{\link{glmm}}.}
}

\value{
 \item{se}{The standard errors for the Monte Carlo maximum likelihood estimates returned from \code{\link{glmm}}}
}

\author{Christina Knudson}

\details{With maximum likelihood performed by Monte Carlo likelihood approximation, there are two sources of variability: there is variability from sample to sample and from Monte Carlo sample (of generated random effects) to Monte Carlo sample. The first source of variability (from sample to sample) is measured using standard error, which appears with the point estimates in the \code{summary} tables. The second source of variability is due to the Monte Carlo randomness, and this is measured by the Monte Carlo standard error. 
}





\seealso{
\code{\link{glmm}} for model fitting.

\code{\link{mcse}} for calculating Monte Carlo standard error.
}
\examples{
library(glmm)
data(BoothHobert)
set.seed(1234)
mod <- glmm(y~0+x1, list(y~0+z1), varcomps.names=c("z1"), 
data=BoothHobert, family.glmm=bernoulli.glmm, m=100, doPQL=TRUE)

se(mod)


}
\keyword{models}
