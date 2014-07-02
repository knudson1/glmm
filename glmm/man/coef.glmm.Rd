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
data(cbpp2)
set.seed(1234)
herd.glmm1<-glmm(Y~period,list(~0+herd), varcomps.names=c("herd"),data=cbpp2, family.glmm=bernoulli.glmm,  m=10^2 )
coef(herd.glmm1)


}

\keyword{generalized linear mixed model}
