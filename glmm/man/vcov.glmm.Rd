\name{vcov.glmm}
\alias{vcov.glmm}

\title{Variance-Covariance Matrix}
\description{A function that calculates the variance-covariance matrix for the Monte Carlo maximum likelihood estimates returned from \code{\link{glmm}}.

}
\usage{
\method{vcov}{glmm}(object,\dots)
}

\arguments{
  \item{object}{An object of class \code{glmm} usually created using \code{\link{glmm}}.}
  \item{\dots}{further arguments passed to or from other methods.}
}

\value{
 \item{vcov }{The variance-covariance matrix for the parameter estimates}
}

\author{Christina Knudson}

\seealso{
\code{\link{glmm}} for model fitting.
}
\examples{
library(glmm)
set.seed(1234)
data(salamander)
#To get more accurate answers for this model, use m=10^4 or 10^5
# and doPQL=TRUE.
m<-10
sal<-glmm(Mate~0+Cross,random=list(~0+Female,~0+Male),varcomps.names=c("F","M"),
 data=salamander,family.glmm=bernoulli.glmm,m=m,debug=TRUE,doPQL=FALSE,cores=2)
vcov(sal)


}

\keyword{generalized linear mixed model}
