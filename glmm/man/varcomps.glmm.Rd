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
data(BoothHobert)
set.seed(1234)
mod <- glmm(y~0+x1, list(y~0+z1), varcomps.names=c("z1"), 
data=BoothHobert, family.glmm=bernoulli.glmm, m=100, doPQL=TRUE)

varcomps(mod)


}

\keyword{generalized linear mixed model}
