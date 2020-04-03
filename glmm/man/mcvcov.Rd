\name{mcvcov}
\alias{mcvcov}

\title{Monte Carlo Variance Covariance Matrix}
\description{A function that calculates the Monte Carlo variance covariance matrix for the Monte Carlo maximum likelihood estimates returned from \code{\link{glmm}}.

}
\usage{mcvcov(object)}

\arguments{
  \item{object}{An object of class \code{glmm} usually created using \code{\link{glmm}}.}
}

\value{
 \item{mcvcov}{The Monte Carlo variance covariance matrix for the Monte Carlo maximum likelihood estimates returned from \code{\link{glmm}}}
}

\author{Christina Knudson}

\details{With maximum likelihood performed by Monte Carlo likelihood approximation, there are two sources of variability: there is variability from sample to sample and from Monte Carlo sample (of generated random effects) to Monte Carlo sample. The first source of variability (from sample to sample) is measured using standard error, which appears with the point estimates in the \code{summary} tables. The second source of variability is due to the Monte Carlo randomness, and this is measured by the Monte Carlo standard error. 

A large entry in Monte Carlo variance covariance matrix indicates the Monte Carlo sample size \code{m} is too small.

The square root of the diagonal elements represents the Monte Carlo standard errors. The off-diagonal entries represent the Monte Carlo covariance.}



\references{
Geyer, C. J. (1994)  On the convergence of Monte Carlo maximum likelihood calculations.  \emph{Journal of the Royal Statistical Society, Series B},  \bold{61}, 261--274.

Knudson, C. (2016)  Monte Carlo likelihood approximation for generalized linear mixed models. PhD thesis, University of Minnesota.  \url{http://hdl.handle.net/11299/178948}}

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
\donttest{mcvcov(sal)}


}

\keyword{generalized linear mixed model}
