\name{mcse}
\alias{mcse}

\title{Monte Carlo Standard Error}
\description{A function that calculates the Monte Carlo standard error for the Monte Carlo maximum likelihood estimates returned from \code{\link{glmm}}.

}
\usage{mcse(object)}

\arguments{
  \item{object}{An object of class \code{glmm} usually created using \code{\link{glmm}}.}
}

\value{
 \item{mcse }{The Monte Carlo standard errors for the Monte Carlo maximum likelihood estimates returned from \code{\link{glmm}}}
}

\author{Christina Knudson}

\details{With maximum likelihood performed by Monte Carlo likelihood approximation, there are two sources of variability: there is variability from sample to sample and from Monte Carlo sample (of generated random effects) to Monte Carlo sample. The first source of variability (from sample to sample) is measured using standard error, which appears with the point estimates in the \code{summary} tables. The second source of variability is due to the Monte Carlo randomness, and this is measured by the Monte Carlo standard error. 

A large Monte Carlo standard error indicates the Monte Carlo sample size \code{m} is too small.}

\references{
Geyer, C. J. (1994)  On the convergence of Monte Carlo maximum likelihood calculations.  \emph{Journal of the Royal Statistical Society, Series B}, \bold{61}, 261--274.

Knudson, C. (2016)  Monte Carlo likelihood approximation for generalized linear mixed models. PhD thesis, University of Minnesota. \url{http://hdl.handle.net/11299/178948}}



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
mcse(sal)


}

\keyword{generalized linear mixed model}
