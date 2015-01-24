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
data(cbpp2)
set.seed(1234)
herd.glmm1<-glmm(Y~period,list(~0+herd), varcomps.names=c("herd"),data=cbpp2, family.glmm=bernoulli.glmm,  m=10^2 )
varcomps(herd.glmm1)


}

\keyword{generalized linear mixed model}
