\name{glmm}
\alias{glmm}
\title{Fitting Generalized Linear Mixed Models using MCML}
\description{This function fits generalized linear mixed models (GLMMs) by approximating the likelihood with ordinary Monte Carlo, then maximizing the approximated likelihood. }
\usage{glmm(fixed, random, varcomps.names, data, family.glmm, m, 
varcomps.equal, weights=NULL, doPQL = TRUE,debug=FALSE, p1=1/3,p2=1/3, p3=1/3,
rmax=1000,iterlim=1000, par.init, zeta=5, cluster=NULL)}


\arguments{
  \item{fixed}{an object of class "\code{formula}" (or one that can be coerced to that class): a symbolic description of the model to be fitted. The details of model specification are given under  "Details."}
  \item{random}{an object of class "\code{formula}" (or one that can be coerced to that class): a symbolic description of the model to be fitted. The details of model specification are given under "Details."}
  \item{varcomps.names}{The names of the distinct variance components in order of \code{varcomps.equal}.}
  \item{data}{an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. If not found in data, the variables are taken from environment(formula), typically the environment from which \code{glmm} is called.}
  \item{family.glmm}{The name of the family. Must be class \code{glmm.family}. Current options are \code{\link{bernoulli.glmm}}, \code{\link{poisson.glmm}}, and \code{\link{binomial.glmm}}.}
  \item{m}{The desired Monte Carlo sample size. See a note in under "Details."}
  \item{varcomps.equal}{An optional vector with elements 1 through the number of distinct variance components. Denotes variance components are to be set equal by assigning them the same integer. The length of varcomps.equal must be equal to the length of the list of random effects formulas. If omitted, varcomps.equal assumes no variance component should be set equal.}
  \item{weights}{An optional vector with length equal to the number of rows of each element of the \code{random} list. This argument makes the specified observations more or less informative for the model.}
  \item{doPQL}{logical. If \code{TRUE}, PQL estimates are used in the importance sampling distribution. If FALSE, the importance sampling distribution will use 0 for the fixed effects and 1 for the variance components. For advanced users, since \code{glmm} is generally more efficient when \code{doPQL=TRUE}.}
  \item{debug}{logical. If \code{TRUE}, extra output useful for testing will be provided. For advanced users.}
\item{p1}{A probability for mixing the random effects generated from three distributions. \code{p1} is the proportion of random effects from the first distribution specified in "Details." For advanced users.}
\item{p2}{A probability for mixing the random effects generated from three distributions. \code{p2} is the proportion of random effects from the second distribution specified in "Details." For advanced users.}
\item{p3}{A probability for mixing the random effects generated from three distributions. \code{p3} is the proportion of random effects from the third distribution specified in "Details." For advanced users.}
\item{rmax}{The maximum allowed trust region radius. This may be set very large. If set small, the algorithm traces a steepest ascent path. This is an argument for \code{trust}.}
\item{iterlim}{A positive integer specifying the maximum number of trust iterations to be performed before the trust program is terminated. This is an argument for \code{trust}.}
\item{par.init}{An optional argument. A single vector that specifies the initial values of the fixed effects and variance components. The parameters should be inputted in the order that \code{\link{summary.glmm}} outputs them, with fixed effects followed by variance components. }
\item{zeta}{A scalar that specifies the degrees of freedom for the t-distribution from which random effects are generated.}
\item{cluster}{An optional argument. A cluster created by the user to enable computations to be done in parallel. See "Details" for more information on the default value.}
}
\details{
Let \eqn{\beta}{beta} be a vector of fixed effects and let \eqn{u}{u} be a vector of random effects. Let \eqn{X}{X} and \eqn{Z}{Z} be design matrices for the fixed and random effects, respectively. The random effects are assumed to be normally distributed with mean 0 and variance matrix \eqn{D}{D}, where \eqn{D}{D} is diagonal with entries from the unknown vector \eqn{\nu}{nu}. Letting \eqn{g}{g} be the link function, \eqn{g(\mu)=X \beta+ZU}{g(mu)=X beta+ ZU}. If the response type is Bernoulli or Binomial, then the logit function is the link; if the response type is Poisson, then the natural logarithm is the link function.

Models for glmm are specified symbolically. A typical fixed effects model has the form \code{response ~ terms} where \code{response} is the (numeric) response vector and \code{terms} is a series of terms which specifies a linear predictor for response. A terms specification of the form \code{first + second} indicates all the terms in \code{first} together with all the terms in \code{second} with duplicates removed. 

A specification of the form \code{first:second} indicates the set of terms obtained by taking the interactions of all terms in \code{first} with all terms in \code{second}. The specification \code{first*second} indicates the cross of \code{first} and \code{second}. This is the same as \code{first + second + first:second}.

The terms in the formula will be re-ordered so that main effects come first, followed by the interactions, all second-order, all third-order and so on: to avoid this, pass a \code{terms} object as the formula. 

If you choose \code{binomial.glmm} as the \code{family.glmm}, then your response should be a two-column matrix: the first column reports the number of successes and the second reports the number of failures.

The random effects for glmm are also specified symbolically. The random effects model specification is typically a list. Each element of the \code{random} list has the form \code{response ~ 0 + term}. The 0 centers the random effects at 0. If you want your random effects to have a nonzero mean, then include that term in the fixed effects.  Each variance component must have its own formula in the list.

To set some variance components equal to one another, use the \code{varcomps.equal} argument. The argument \code{varcomps.equal} should be a vector whose length is equal to the length of the random effects list. The vector should contain positive integers, and the first element of the \code{varcomps.equal} should be 1. To set variance components equal to one another, assign the same integer to the corresponding elements of \code{varcomps.equal}. For example, to set the first and second variance components equal to each other, the first two elements of \code{varcomps.equal} should be 1. If \code{varcomps.equal} is omitted, then the variance components are assumed to be distinct. 

Each distinct variance component should have a name. The length of \code{varcomps.names} should be equal to the number of distinct variance components. If \code{varcomps.equal} is omitted, then the length of \code{varcomps.names} should be equal to the length of \code{random}.

Monte Carlo likelihood approximation relies on an importance sampling distribution. Though infinitely many importance sampling distributions should yield the correct MCMLEs eventually, the importance sampling distribution used in this package was chosen to reduce the computation cost. When \code{doPQL} is \code{TRUE}, the importance sampling distribution relies on PQL estimates (as calculated in this package). When \code{doPQL} is \code{FALSE}, the random effect estimates in the distribution are taken to be 0, the fixed effect estimates are taken to be 0, and the variance component estimates are taken to be 1.

This package's importance sampling distribution is a mixture of three distributions: a t centered at 0 with scale matrix determined by the PQL estimates of the variance components and with \code{zeta} degrees of freedom, a normal distribution centered at the PQL estimates of the random effects and with a variance matrix containing the PQL estimates of the variance components, and a normal distribution centered at the PQL estimates of the random effects and with a variance matrix based on the Hessian of the penalized log likelihood. The first component is included to guarantee the gradient of the MCLA has a central limit theorem. The second component is included to mirror our best guess of the distribution of the random effects. The third component is included so that the numerator and the denominator are similar when calculating the MCLA value.

The Monte Carlo sample size \code{m} should be chosen as large as possible. You may want to run the model a couple times to begin to understand the variability inherent to Monte Carlo. There are no hard and fast rules for choosing \code{m}, and more research is needed on this area. For a general idea, I believe the \code{BoothHobert} model produces stable enough estimates at \eqn{m=10^3} and the \code{salamander} model produces stable enough estimates at \eqn{m=10^5}, as long as \code{doPQL} is \code{TRUE}.

To decrease the computation time involved with larger Monte Carlo sample sizes, a parallel computing component has been added to this package. Generally, the larger the Monte Carlo sample size \code{m}, the greater the benefit of utilizing additional cores. By default, \code{glmm} will create a cluster that uses a single core. This forces all computations to be done sequentially rather than simultaneously. 

To see the summary of the model, use summary().


}

\value{

\code{glmm} returns an object of class \code{glmm} is a list containing at least the following components: 
\item{beta }{A vector of the Monte Carlo maximum likelihood estimates (MCMLEs) for the fixed effects.}
 \item{nu }{A vector of the Monte Carlo maximum likelihood estimates for the variance components.}
 \item{loglike.value}{The Monte Carlo log likelihood evaluated at the MCMLEs \code{beta} and \code{nu}.}
 \item{loglike.gradient}{The Monte Carlo log likelihood gradient vector at the MCMLEs \code{beta} and \code{nu}.}
 \item{loglike.hessian}{The Monte Carlo log likelihood Hessian matrix at the MCMLEs \code{beta} and \code{nu}.}
\item{mod.mcml}{A list containing the fixed effect design matrix, the list of random effect design matrices, the response. and the number of trials (for the Binomial family).}
\item{call}{The call.}
\item{fixedcall}{The fixed effects call.}
\item{randcall}{The random effects call.}
\item{x}{The design matrix for the fixed effects.}
\item{y}{The response vector.}
\item{z}{The design matrix for the random effects.}
\item{family.glmm}{The name of the family. Must be class \code{glmm.family}.}
\item{varcomps.names}{The vector of names for the distinct variance components.}
\item{varcomps.equal}{The vector denoting equal variance components.}
\item{umat}{A matrix with \code{m} rows. Each row is a vector of random effects generated from the importance sampling distribution.}
\item{pvec}{A vector containing \code{p1}, \code{p2}, and \code{p3}.}
\item{beta.pql}{PQL estimate of \eqn{\beta}{beta}, when \code{doPQL} is \code{TRUE}.}
\item{nu.pql}{PQL estimate of \eqn{\nu}{nu}, when \code{doPQL} is \code{TRUE}.}
\item{u.pql}{PQL predictions of the random effects.}
\item{zeta}{The number of degrees of freedom used in the t component of the importance sampling distribution.}
\item{cluster}{A cluster created by the user for use during the approximation of the log-likelihood value, gradient and hessian.}
  \item{debug}{If \code{TRUE} extra output useful for testing.}
  The function \code{\link{summary}} (i.e., \code{\link{summary.glmm}}) can
  be used to obtain or print a summary of the results. The generic accessor function \code{\link{coef}} (i.e., \code{\link{coef.glmm}}) can be used to extract the coefficients.
}
\references{
Geyer, C. J. (1994)  On the convergence of Monte Carlo maximum likelihood calculations.  \emph{Journal of the Royal Statistical Society, Series B}, \bold{61}, 261--274.

Geyer, C. J. and Thompson, E. (1992)  Constrained Monte Carlo maximum likelihood for dependent data.  \emph{Journal of the Royal Statistical Society, Series B}, \bold{54}, 657--699.

Knudson, C. (2016). Monte Carlo likelihood approximation for generalized linear mixed models. PhD thesis, University of Minnesota.  \url{http://hdl.handle.net/11299/178948}

Sung, Y. J. and Geyer, C. J. (2007)  Monte Carlo likelihood inference for missing data models.  \emph{Annals of Statistics}, \bold{35}, 990--1011.}

\author{Christina Knudson}
%\note{}





\examples{
#First, using the basic Booth and Hobert dataset
#to fit a glmm with a logistic link, one variance component,
#one fixed effect, and an intercept of 0. The Monte Carlo
#sample size is 100 to save time.
library(glmm)
data(BoothHobert)
set.seed(1234)
\donttest{mod.mcml1 <- glmm(y~0+x1, list(y~0+z1), varcomps.names=c("z1"), data=BoothHobert,}
\donttest{family.glmm=bernoulli.glmm, m=100, doPQL=TRUE, cluster=clust)}
\donttest{mod.mcml1$beta}
\donttest{mod.mcml1$nu}
\donttest{summary(mod.mcml1)}
\donttest{coef(mod.mcml1)}

#Next, a model setting two variance components equal.
data(Booth2)
set.seed(1234)
\donttest{mod.mcml3<-glmm(y~0+x1, list(y~0+z1,~0+z2), varcomps.names=c("z"), }
\donttest{varcomps.equal=c(1,1), data=Booth2, family.glmm=bernoulli.glmm, m=100)}
\donttest{summary(mod.mcml3)}

#Now, a model with crossed random effects. There are two distinct
#variance components. To get more accurate answers for this model,
#use a larger Monte Carlo sample size, such as m=10^4 or 10^5
#and doPQL=TRUE.
set.seed(1234)
data(salamander)
m<-100
\donttest{sal<-glmm(Mate~0+Cross,random=list(~0+Female,~0+Male),varcomps.names=c("F","M"),}
\donttest{data=salamander,family.glmm=bernoulli.glmm,m=m,debug=TRUE,cluster=clust)}
\donttest{summary(sal)}

#The above model (sal) can be redone with binomial.glmm
set.seed(1234)
\donttest{sal<-glmm(cbind(Mate,1-Mate)~0+Cross,random=list(~0+Female,~0+Male),varcomps.names=c("F","M"),}
\donttest{data=salamander,family.glmm=binomial.glmm,m=m,debug=TRUE,cluster=clust)}

}

\keyword{generalized linear mixed model }
\keyword{Monte Carlo }
\keyword{maximum likelihood}
\keyword{likelihood approximation}
