\name{logLik.glmm}
\alias{logLik.glmm}

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
set.seed(1234)
data(salamander)

m<-1000
\donttest{sal<-glmm(Mate~0+Cross,random=list(~0+Female,~0+Male),varcomps.names=c("F","M"),}
\donttest{data=salamander,family.glmm=bernoulli.glmm,m=m,debug=TRUE,cores=2)}
\donttest{logLik(sal)}


}

\keyword{generalized linear mixed model}
