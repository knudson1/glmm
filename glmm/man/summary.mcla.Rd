\name{summary.mcla}
\alias{summary.mcla}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
summary.mcla(mod.mcml)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mod.mcml}{
%%     ~~Describe \code{mod.mcml} here~~
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
function (mod.mcml) 
{
    stopifnot(inherits(mod.mcml, "mcla"))
    fixedcall <- mod.mcml$fixedcall
    randcall <- mod.mcml$randcall
    call <- mod.mcml$call
    x <- mod.mcml$x
    y <- mod.mcml$y
    z <- mod.mcml$z
    beta <- mod.mcml$beta
    nbeta <- length(beta)
    hessian <- mod.mcml$likelihood.hessian
    if (det(hessian) == 0) {
        warning(paste("estimated Fisher information matrix not positive", 
            "definite, making all standard errors infinite"))
        all.ses <- rep(Inf, nrow(hessian))
    }
    else {
        all.ses <- sqrt(diag(solve(-hessian)))
    }
    beta.se <- all.ses[1:nbeta]
    zval <- beta/beta.se
    coefmat <- cbind(beta, beta.se, zval, 2 * pnorm(abs(zval), 
        lower.tail = F))
    colnames(coefmat) <- c("Estimate", "Std. Error", "z value", 
        "Pr(>|z|)")
    rownames(coefmat) <- colnames(mod.mcml$x)
    nu <- mod.mcml$nu
    nu.se <- all.ses[-(1:nbeta)]
    nuzval <- nu/nu.se
    nucoefmat <- cbind(nu, nu.se, nuzval, pnorm(abs(nuzval), 
        lower.tail = F))
    colnames(nucoefmat) <- c("Estimate", "Std. Error", "z value", 
        "Pr(>|z|)/2")
    rownames(nucoefmat) <- mod.mcml$varcomps.names
    return(structure(list(x = x, y = y, z = z, coefmat = coefmat, 
        fixedcall = fixedcall, randcall = randcall, family.mcml = mod.mcml$family.mcml, 
        call = call, nucoefmat = nucoefmat), class = "summary.mcla"))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
