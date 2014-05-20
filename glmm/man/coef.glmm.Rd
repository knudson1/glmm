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
herd.glmm1<-glmm(Y~period,list(~0+herd), family=bernoulli.glmm, data=cbpp2, m=10^2, varcomps.names=c("herd"))
coef(herd.glmm1)


}

\keyword{generalized linear mixed model}
