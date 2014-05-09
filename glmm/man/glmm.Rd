\name{glmm}
\alias{glmm}
\alias{mcml}
\title{Fitting Generalized Linear Mixed Models using MCML}
\description{This function fits generalized linear mixed models (GLMMs) by approximating the likelihood with ordinary Monte Carlo, then maximizing the approximated likelihood. }
\usage{glmm(fixed, random, varcomps.names, data, family.mcml, m, 
varcomps.equal, doPQL = TRUE)}


\arguments{
  \item{fixed}{an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted. The details of model specification are given under  "Details."}
  \item{random}{an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted. The details of model specification are given under Details.}
  \item{varcomps.names}{The names of the distinct variance components in order of varcomps.equal.}
  \item{data}{an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. If not found in data, the variables are taken from environment(formula), typically the environment from which lm is called.}
  \item{family.mcml}{The name of the family. Must be class mcml.family.}
  \item{m}{The desired Monte Carlo sample size.}
  \item{varcomps.equal}{An optional vector with elements 1 through the number of distinct variance components. Denotes variance components are to be set equal by assigning them the same integer. The length of varcomps.equal must be equal to the length of the list of random effects formulas. If omitted, varcomps.equal assumes no variance component should be set equal.}
  \item{doPQL}{logical. If TRUE, PQL estimates are used for the importance sampling distribution. If FALSE, the importance sampling distribution is standard normal.}
}
\details{
Models for glmm are specified symbolically. A typical fixed effects model has the form \code{response ~ terms} where \code{response} is the (numeric) response vector and \code{terms} is a series of terms which specifies a linear predictor for response. A terms specification of the form \code{first + second} indicates all the terms in \code{first} together with all the terms in \code{second} with duplicates removed. 

A specification of the form \code{first:second} indicates the set of terms obtained by taking the interactions of all terms in \code{first} with all terms in \code{second}. The specification \code{first*second} indicates the cross of \code{first} and \code{second}. This is the same as \code{first + second + first:second}.

The terms in the formula will be re-ordered so that main effects come first, followed by the interactions, all second-order, all third-order and so on: to avoid this, pass a \code{terms} object as the formula. 

The random effects for glmm are also specified symbolically. The random effects model specification is typically a list. Each element of the \code{random} list has the form \code{response ~ 0 + term}. The 0 is needed to center the random effects at 0. If you want your random effects to have a nonzero mean, then include that term in the fixed effects.  Each variance component must have its own formula in the list.

To set some variance components equal to one another, use the \code{varcomps.equal} argument. The argument \code{varcomps.equal} should be a vector whose length is equal to the length of the random effects list. The vector should contain nonnegative integers, and the first element of the \code{varcomps.equal} should be 1. To set variance components equal to one another, assign the same integer to the corresponding elements of \code{varcomps.equal}. For example, to set the first and second variance components equal to each other, the first two elements of \code{varcomps.equal} should be 1. If \code{varcomps.equal} is omitted, then the variance components are assumed to be distinct. 

Each distinct variance component should have a name. The length of \code{varcomps.names} should be equal to the number of distinct variance components. If \code{varcomps.equal} is omitted, then the length of \code{varcomps.names} should be equal to the length of \code{random}.

To see the summary of the model use summary().
}

\value{

\code{glmm} returns an object of class "mcla" is a list containing at least the following components: 
\item{beta }{A vector of the Monte Carlo maximum likelihood estimates for the fixed effects.}
 \item{nu }{A vector of the Monte Carlo maximum likelihood estimates for the variance components.}
 \item{likelihood.value}{The Monte Carlo maximum likelihood evaluated at the MCMLEs beta and nu.}
 \item{likelihood.gradient}{The Monte Carlo maximum likelihood vector at the MCMLEs beta and nu.}
 \item{likelihood.hessian}{The Monte Carlo maximum likelihood hessian matrix at the MCMLEs beta and nu.}
 \item{trust.converged}{if TRUE the final iterate was deemed optimal by the termination criteria specified by trust. The MCMLE was identified.}
\item{beta.pql }{A vector of the beta values used in the importance sampling distribution.}
 \item{nu.pql }{A vector of the nu values used in the importance sampling distribution.}
\item{mod.mcml}{A list containing the fixed effect design matrix, the list of random effect design matrices, and the response.}
\item{trust.argpath}{The sequence of iterates in the trust optimization, not including the final iterate.}
\item{fixedcall}{The fixed effects call.}
\item{fixedcall}{The random effects call.}
\item{x}{The design matrix for the fixed effects.}
\item{y}{The response vector.}
\item{z}{The design matrix for the random effects.}
\item{family.mcml}{The name of the family. Must be class mcml.family.}
\item{call}{The call}
\item{mat}{The matrix of random effects generated according to the importance sampling distribution with u.pql and nu.pql as parameters.}
\item{varcomps.names}{The vector of names for the distinct variance components.}
\item{varcomps.equal}{The vector denoting equal variance components.}
\item{u.pql}{The vector of random effects as optimized by PQL. This vector serves as the mean vector for the importance sampling distribution.}
  The function \code{\link{summary}} (i.e., \code{\link{summary.mcla}}) can
  be used to obtain or print a summary of the results. The generic accessor function \code{\link{coef}} (i.e., \code{\link{coef.mcla}}) can be used to extract the coefficients.
}
\references{
Geyer, C. J. (1994). On the convergence of Monte Carlo maximum likelihood calculations. Journal of the Royal Statistical Society, Series B, 61:261--274.

Geyer, C. J. and Thompson, E. (1992). Constrained Monte Carlo maximum likelihood for dependent data. Journal of the Royal Statistical Society, Series B, 54:657--699.

Sung, Y. J. and Geyer, C. J. (2007). Monte Carlo likelihood inference for missing data models. Annals of Statistics, 35:990--1011.}

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
mod.mcml1<-glmm(y~0+x1,list(y~0+z1),varcomps.names=c("z1"),data=BoothHobert,
family.mcml=bernoulli.mcml,m=100,doPQL=TRUE)
summary(mod.mcml1)
coef(mod.mcml1)

#Next, a similar model that has two distinct variance components.
data(Booth2)
set.seed(1234)
mod.mcml2<-glmm(y~0+x1,list(y~0+z1,~0+z2),varcomps.names=c("z1","z2"),
data=Booth2,family.mcml=bernoulli.mcml,m=100,doPQL=TRUE)
summary(mod.mcml2)

#Next, a similar model but we set the variance components equal.
data(Booth2)
set.seed(1234)
mod.mcml3<-glmm(y~0+x1,list(y~0+z1,~0+z2),varcomps.names=c("z"), 
varcomps.equal=c(1,1), data=Booth2,family.mcml=bernoulli.mcml,
m=100,doPQL=FALSE)
summary(mod.mcml3)

}

\keyword{generalized linear mixed model }
\keyword{Monte Carlo }
\keyword{maximum likelihood}
\keyword{likelihood approximation}
