\name{varcomps}
\alias{varcomps}

\title{Extract Model Variance Components}
\description{A function that extracts the variance components returned from \code{\link{glmm}}.

}
\usage{varcomps(object,\dots)}


\arguments{
  \item{object}{An object of class \code{glmm} usually created using \code{\link{glmm}}.}
  \item{\dots}{further arguments passed to or from other methods.}
}

\value{
 \item{varcomps }{A vector of variance component estimates }
}

\author{Christina Knudson}

\seealso{
\code{\link{glmm}} for model fitting.
\code{\link{coef.glmm}} for fixed effects coefficients.
}
\examples{
library(glmm)
set.seed(1234)
data(salamander)
#To get more accurate answers for this model, use m=10^4 or 10^5
#and doPQL=TRUE.
m<-10
sal<-glmm(Mate~0+Cross,random=list(~0+Female,~0+Male),varcomps.names=c("F","M"),
data=salamander,family.glmm=bernoulli.glmm,m=m,debug=TRUE,doPQL=FALSE)
varcomps(sal)


}

\keyword{generalized linear mixed model}
