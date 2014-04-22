\name{glmm}
\alias{mcml}
\title{Fitting Generalized Linear Mixed Models using MCML}
\description{This function fits generalized linear mixed models (GLMMs) by approximating the likelihood with ordinary Monte Carlo, then maximizing the approximated likelihood. }
\usage{glmm(fixed, random, varcomps.names, data, family.mcml, m, varcomps.equal, doPQL = T)}
\usage{summary(mod.mcml)}

\arguments{
  \item{fixed}{an object of class \dquote formula \dquote (or one that can be coerced to that class): a symbolic description of the model to be fitted. The details of model specification are given under \squote Details. \squote}
  \item{random}{an object of class \dquote formula \dquote (or one that can be coerced to that class): a symbolic description of the model to be fitted. The details of model specification are given under DetailsÕ.}
  \item{varcomps.names}{The names of the distinct variance components in order of varcomps.equal.}
  \item{data}{an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. If not found in data, the variables are taken from environment(formula), typically the environment from which lm is called.}
  \item{family.mcml}{The name of the family. Must be class mcml.family.}
  \item{m}{The desired Monte Carlo sample size.}
  \item{varcomps.equal}{An optional vector with elements 1 through the number of distinct variance components. Denotes variance components are to be set equal by assigning them the same integer. The length of varcomps.equal must be equal to the length of the list of random effects formulas. If omitted, varcomps.equal assumes no variance component should be set equal.}
  \item{doPQL}{logical. If TRUE, PQL estimates are used for the importance sampling distribution. If FALSE, the importance sampling distribution is standard normal.}
}
\details{
Models for glmm are specified symbolically. A typical fixed effects model has the form response ~ terms where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response. A terms specification of the form first + second indicates all the terms in first together with all the terms in second with duplicates removed. A specification of the form first:second indicates the set of terms obtained by taking the interactions of all terms in first with all terms in second. The specification first*second indicates the cross of first and second. This is the same as first + second + first:second.

The random effects for glmm are also specified symbolically. The random effects model specification is typically a list. Each element of the list has the form response ~ 0 + term. The 0 is needed to center the random effects at 0. If you want your random effects to have a nonzero mean, then include that term in the fixed effects.

If you want to set some of your variance components equal to one another, use the varcomps.equal argument. The argument varcomps.equal should be a vector whose length is equal to the length of the random effects list. If varcomps.equal is omitted, then the variance components are assumed to be distinct. 

To see the summary of the model use summary().
}

\value{
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

}
\references{
Geyer, C. J. (1994). On the convergence of Monte Carlo maximum likelihood calculations. Journal of the Royal Statistical Society, Series B, 61:261–274.

Geyer, C. J. and Thompson, E. (1992). Constrained Monte Carlo maximum likelihood for dependent data. Journal of the Royal Statistical Society, Series B, 54:657–699.

Sung, Y. J. and Geyer, C. J. (2007). Monte Carlo likelihood inference for missing data models. Annals of Statistics, 35:990–1011.}

\author{Christina Knudson}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}

\examples{
library(glmm)
data(BoothHobert)
mod.mcml<-glmm(y~0+x1,list(y~0+z1),varcomps.names=c("z1"),data=dat,family.mcml=bernoulli.mcml,m=100,doPQL=T)
summary(mod.mcml)

data(Booth2)

}

\keyword{generalized linear mixed model }
\keyword{Monte Carlo }
\keyword{maximum likelihood}
\keyword{likelihood approximation}
