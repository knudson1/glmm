\name{coef.glmm}
\alias{coef.glmm}
\concept{generalized linear mixed model}

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
m<-1000
\donttest{sal<-glmm(Mate~0+Cross,random=list(~0+Female,~0+Male),varcomps.names=c("F","M"),}
\donttest{data=salamander,family.glmm=bernoulli.glmm,m=m,debug=TRUE,doPQL=FALSE)}
\donttest{coef(sal)}

}

\keyword{regression}
\keyword{models}

