\name{mcml}
\alias{mcml}
\title{Fitting Generalized Linear Mixed Models using MCML}
\description{This function fits generalized linear mixed models (GLMMs) by approximating the likelihood with ordinary Monte Carlo, then maximizing the approximated likelihood. }
\usage{mcml(fixed, random, varcomps.names, data, family.mcml, m, varcomps.equal, doPQL = T)}
\usage{summary(mod.mcml)}
\arguments{
  \item{fixed}{an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted. The details of model specification are given under ÔDetailsÕ.}
  \item{random}{an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted. The details of model specification are given under ÔDetailsÕ.}
  \item{varcomps.names}{The names of the distinct variance components in order of varcomps.equal.}
  \item{data}{an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. If not found in data, the variables are taken from environment(formula), typically the environment from which lm is called.}
  \item{family.mcml}{The name of the family. Must be class mcml.family.}
  \item{m}{The desired Monte Carlo sample size.}
  \item{varcomps.equal}{An optional vector with elements 1 through the number of distinct variance components. Denotes variance components are to be set equal by assigning them the same integer. The length of varcomps.equal must be equal to the length of the list of random effects formulas. If omitted, varcomps.equal assumes no variance component should be set equal.}
  \item{doPQL}{logical. If TRUE, PQL estimates are used for the importance sampling distribution. If FALSE, the importance sampling distribution is standard normal.}
}
\details{
Models for mcml are specified symbolically. A typical model has the form response ~ terms where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response. A terms specification of the form first + second indicates all the terms in first together with all the terms in second with duplicates removed. A specification of the form first:second indicates the set of terms obtained by taking the interactions of all terms in first with all terms in second. The specification first*second indicates the cross of first and second. This is the same as first + second + first:second.

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
%% ~put references to the literature/web site here ~
}
\author{Christina Knudson}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (fixed, random, varcomps.names, data, family.mcml, m, 
    varcomps.equal, doPQL = T) 
{
    if (missing(varcomps.names)) 
        stop("Names for the variance components must be supplied through varcomps.names")
    if (is.vector(varcomps.names) != 1) 
        stop("varcomps.names must be a vector")
    if (missing(varcomps.equal)) {
        varcomps.equal <- c(1:length(varcomps.names))
    }
    call <- match.call()
    stopifnot(inherits(fixed, "formula"))
    if (missing(data)) {
        barf <- lm(fixed, method = "model.frame")
    }
    else {
        stopifnot(inherits(data, "data.frame"))
        barf <- lm(fixed, data = data, method = "model.frame")
    }
    x <- model.matrix(fixed, data = barf)
    y <- model.response(barf)
    randcall <- random
    if (!is.list(random)) 
        random <- list(random)
    for (irandom in seq(along = random)) {
        r <- random[[irandom]]
        stopifnot(inherits(r, "formula"))
        if (missing(data)) {
            barf2 <- lm(r, method = "model.frame")
        }
        else {
            stopifnot(inherits(data, "data.frame"))
            barf2 <- lm(r, data = data, method = "model.frame")
        }
        random[[irandom]] <- model.matrix(r, data = barf2)
        if (length(y) != nrow(random[[irandom]])) 
            stop("Fixed and random effect model matrices should have same number of rows")
    }
    if (is.numeric(varcomps.equal) == F) 
        stop("varcomps.equal must be a vector containing numbers to indicate which variance components are equal.")
    if (length(varcomps.equal) != length(random)) {
        stop("The length of varcomps.equal must be equal to the length of the random-effects call.")
    }
    if (length(unique(varcomps.equal)) != length(varcomps.names)) {
        stop("You must name each unique variance component. Check varcomps.names and varcomps.equal.")
    }
    if (min(varcomps.equal) != 1) 
        stop("The vector varcomps.equal must contain numbers starting at 1 to denote which variance components are equal.")
    levs <- ordered(unique(varcomps.equal))
    z <- list()
    for (i in 1:length(levs)) {
        if (levs[i] != i) 
            stop("The numbers in the vector varcomps.equal must be consecutive. You must start at 1 and then each entry must be the next consecutive number or a repeat of a previous number.")
        these <- varcomps.equal == i
        thesemats <- random[these]
        z[[i]] <- do.call(cbind, thesemats)
    }
    names(z) <- varcomps.names
    mod.mcml <- structure(list(x = x, z = z, y = y), class = "bar")
    mod.mcml
    if (doPQL == T) {
        pql.out <- pql(mod.mcml, family.mcml)
        s.pql <- pql.out$s
        sigma.pql <- pql.out$sigma
        nu.pql <- sigma.pql^2
        beta.pql <- pql.out$beta
    }
    if (doPQL == F) {
        nrand <- lapply(mod.mcml$z, ncol)
        nrandom <- unlist(nrand)
        totnrandom <- sum(nrandom)
        s.pql <- rep(0, totnrandom)
        nu.pql <- rep(1, length(mod.mcml$z))
        beta.pql <- rep(1, ncol(mod.mcml$x))
    }
    par.init <- c(pql.out$beta, nu.pql)
    genData <- genRand(sigma.pql, s.pql, mod.mcml$z, m)
    umat <- genData$u
    u.star <- genData$u.star
    trust.out <- trust(objfun, parinit = par.init, rinit = 10, 
        rmax = 10000, iterlim = 100, minimize = F, nbeta = length(beta.pql), 
        nu.pql = nu.pql, umat = umat, mod.mcml = mod.mcml, family.mcml = family.mcml, 
        m = m, u.star = u.star, blather = T)
    beta.trust <- trust.out$argument[1:length(beta.pql)]
    nu.trust <- trust.out$argument[-(1:length(beta.pql))]
    trust.argpath <- trust.out$argpath
    names(beta.trust) <- colnames(mod.mcml$x)
    names(nu.trust) <- varcomps.names
    return(structure(list(beta = beta.trust, nu = nu.trust, likelihood.value = trust.out$value, 
        likelihood.gradient = trust.out$gradient, likelihood.hessian = trust.out$hessian, 
        trust.converged = trust.out$converged, beta.pql = beta.pql, 
        nu.pql = nu.pql, mod.mcml = mod.mcml, trust.argpath = trust.argpath, 
        fixedcall = fixed, randcall = randcall, x = x, y = y, 
        z = random, family.mcml = family.mcml, call = call, umat = umat, 
        varcomps.names = varcomps.names, varcomps.equal = varcomps.equal, 
        u.pql = u.star), class = "mcla"))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
