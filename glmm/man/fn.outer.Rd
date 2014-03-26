\name{fn.outer}
\alias{fn.outer}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
fn.outer(par, beta, s, Y, X, Z, eek, family.mcml)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{par}{
%%     ~~Describe \code{par} here~~
}
  \item{beta}{
%%     ~~Describe \code{beta} here~~
}
  \item{s}{
%%     ~~Describe \code{s} here~~
}
  \item{Y}{
%%     ~~Describe \code{Y} here~~
}
  \item{X}{
%%     ~~Describe \code{X} here~~
}
  \item{Z}{
%%     ~~Describe \code{Z} here~~
}
  \item{eek}{
%%     ~~Describe \code{eek} here~~
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
function (par, beta, s, Y, X, Z, eek, family.mcml) 
{
    sigma <- par
    Aks <- Map("*", eek, sigma)
    A <- addVecs(Aks)
    A <- diag(A)
    nbeta <- length(beta)
    inner.optim <- trust(fn.inner.trust, parinit = c(beta, s), 
        rinit = 5, rmax = 10000, minimize = F, Y = Y, X = X, 
        Z = Z, A = A, nbeta = nbeta, family.mcml = family.mcml)
    parms <- inner.optim$argument
    beta.twid <<- parms[1:nbeta]
    s.twid <<- parms[-(1:nbeta)]
    eta.twid <- (X \%*\% beta.twid + Z \%*\% A \%*\% s.twid)
    family.mcml <- getFamily(family.mcml)
    piece1 <- el(Y, X, eta.twid, family.mcml)$value
    dubya <- family.mcml$cpp(eta.twid)
    W <- diag(as.vector(dubya))
    thatthing <- A \%*\% t(Z) \%*\% W \%*\% Z \%*\% A + diag(nrow(A))
    L <- chol(thatthing)
    piece2a <- sum(log(diag(L)))
    piece2b <- 0.5 * s \%*\% s
    minvalue <- piece2a + piece2b - piece1
    as.vector(minvalue)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
