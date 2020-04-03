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
m<-1000
\donttest{sal<-glmm(Mate~0+Cross,random=list(~0+Female,~0+Male),varcomps.names=c("F","M"),}
\donttest{data=salamander,family.glmm=bernoulli.glmm,m=m,debug=TRUE)}
\donttest{vcov(sal)}


}

\keyword{generalized linear mixed model}
