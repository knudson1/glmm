\name{mcml}
\alias{mcml}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
mcml(fixed, random, varcomps.names, data, family.mcml, m, varcomps.equal, doPQL = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fixed}{
%%     ~~Describe \code{fixed} here~~
}
  \item{random}{
%%     ~~Describe \code{random} here~~
}
  \item{varcomps.names}{
%%     ~~Describe \code{varcomps.names} here~~
}
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{family.mcml}{
%%     ~~Describe \code{family.mcml} here~~
}
  \item{m}{
%%     ~~Describe \code{m} here~~
}
  \item{varcomps.equal}{
%%     ~~Describe \code{varcomps.equal} here~~
}
  \item{doPQL}{
%%     ~~Describe \code{doPQL} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
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
