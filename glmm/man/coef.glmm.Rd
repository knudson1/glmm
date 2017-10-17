\name{coef.glmm}
\alias{coef.glmm}

\title{Extract Model Coefficients}
\description{A function that extracts the fixed effect coefficients returned from \code{\link{glmm}}.

}
\usage{
\method{coef}{glmm}(object,\dots)
}

\arguments{
  \item{object}{An object of class \code{glmm} usually created using \code{\link{glmm}}.}
  \item{\dots}{further arguments passed to or from other methods.}
}

\value{
 \item{coefficients }{A vector of coefficients (fixed effects only)}
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
coef(sal)


}

\keyword{generalized linear mixed model}
