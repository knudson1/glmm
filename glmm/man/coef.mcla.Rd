\name{coef.mcla}
\alias{coef.mcla}
\alias{coefficients.mcla}
\title{Extract Model Coefficients}
\description{A function that extracts the fixed effect coefficients returned from \code{\link{glmm}}.

}
\usage{
coef(mod.mcml)
coefficients(mod.mcml)
}

\arguments{
  \item{mod.mcml}{An object of class "mcla" usually created using \code{\link{glmm}}.
}
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
data(BoothHobert)
mod.mcml1<-glmm(y~0+x1,list(y~0+z1),varcomps.names=c("z1"),
data=BoothHobert,family.mcml=bernoulli.mcml,m=100,doPQL=TRUE)
coef(mod.mcml1)
coefficients(mod.mcml1)

}

\keyword{generalized linear mixed model}
