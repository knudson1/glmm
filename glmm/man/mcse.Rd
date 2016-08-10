\name{mcse}
\alias{mcse}

\title{Monte Carlo Standard Error}
\description{A function that calculates the Monte Carlo standard error for the Monte Carlo maximum likelihood estimates returned from \code{\link{glmm}}.

}
\usage{mcse(object)}

\arguments{
  \item{object}{An object of class \code{glmm} usually created using \code{\link{glmm}}.}
}

\value{
 \item{mcse }{The Monte Carlo standard errors for the Monte Carlo maximum likelihood estimates returned from \code{\link{glmm}}}
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
 data=salamander,family.glmm=bernoulli.glmm,m=m,debug=TRUE,doPQL=FALSE)
mcse(sal)


}

\keyword{generalized linear mixed model}
