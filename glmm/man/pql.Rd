\name{pql}
\alias{pql}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
pql(mod.mcml, family.mcml)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mod.mcml}{
%%     ~~Describe \code{mod.mcml} here~~
}
  \item{family.mcml}{
%%     ~~Describe \code{family.mcml} here~~
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
function (mod.mcml, family.mcml) 
{
    eek <- getEk(mod.mcml$z)
    sigma <- rep(1, length(mod.mcml$z))
    beta <- rep(0, ncol(mod.mcml$x))
    nrand <- lapply(mod.mcml$z, ncol)
    nrandom <- unlist(nrand)
    totnrandom <- sum(nrandom)
    s <- rep(1, totnrandom)
    Y = mod.mcml$y
    X = mod.mcml$x
    Z = do.call(cbind, mod.mcml$z)
    outer.optim <- suppressWarnings(optim(par = sigma, fn = fn.outer, 
        beta = beta, s = s, Y = Y, X = X, Z = Z, eek = eek, family.mcml = family.mcml))
    list(sigma = outer.optim$par, beta = beta.twid, s = s.twid)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
