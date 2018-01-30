\name{confint.glmm}
\alias{confint.glmm}

\title{Calculates Asymptotic Confidence Intervals}
\description{A function that calculates asymptotic confidence intervals for one or more parameters in a model fitted by by \code{\link{glmm}}. Confidence intervals can be calculated for fixed effect parameters and variance components using models.

}
\usage{
\method{confint}{glmm}(object, parm, level, \dots)
}

\arguments{
  \item{object}{An object of class \code{glmm} usually created using \code{\link{glmm}}.}
  \item{parm}{A specification of which parameters are to be given confidence intervals, either a vector of numbers or a vector of names. If missing, all parameters are considered.}
  \item{level}{The confidence level required.}
  \item{\dots}{Additional arguments passed to or from other methods.}
}

\value{A matrix (or vector) with columns giving lower and upper confidence limits for each parameter. These will be labeled as (1-level)/2 and 1-(1-level)/2 in percent. By default, 2.5% and 97.5%.}
}

\author{Christina Knudson}

\seealso{
\code{\link{glmm}} for model fitting.
}
\examples{
library(glmm)
data(BoothHobert)
set.seed(123)
mod.mcml1<-glmm(y~0+x1,list(y~0+z1),varcomps.names=c("z1"), 
data=BoothHobert,family.glmm=bernoulli.glmm,m=10,doPQL=TRUE)
confint(mod.mcml1)

}

\keyword{generalized linear mixed model}
