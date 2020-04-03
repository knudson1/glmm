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
m<-1000
\donttest{sal<-glmm(Mate~0+Cross,random=list(~0+Female,~0+Male),varcomps.names=c("F","M"),}
\donttest{data=salamander,family.glmm=bernoulli.glmm,m=m,debug=TRUE)}
\donttest{varcomps(sal)}


}

\keyword{generalized linear mixed model}
