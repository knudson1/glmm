\name{objfun}
\alias{objfun}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
objfun(par, nbeta, nu.pql, umat, u.star = u.star, mod.mcml, family.mcml)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{par}{
%%     ~~Describe \code{par} here~~
}
  \item{nbeta}{
%%     ~~Describe \code{nbeta} here~~
}
  \item{nu.pql}{
%%     ~~Describe \code{nu.pql} here~~
}
  \item{umat}{
%%     ~~Describe \code{umat} here~~
}
  \item{u.star}{
%%     ~~Describe \code{u.star} here~~
}
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
function (par, nbeta, nu.pql, umat, u.star = u.star, mod.mcml, 
    family.mcml) 
{
    beta <- par[1:nbeta]
    nu <- par[-(1:nbeta)]
    m <- nrow(umat)
    if (sum(nu <= 0) > 0) {
        out <- list(value = -Inf, gradient = rep(1, length(par)), 
            hessian = as.matrix(c(rep(1, length(par)^2)), nrow = length(par)))
        return(out)
    }
    Z = do.call(cbind, mod.mcml$z)
    eta <- b <- rep(0, m)
    lfu <- lfu.twid <- lfyu <- list(rep(c(0, 0, 0), m))
    for (k in 1:m) {
        Uk <- umat[k, ]
        eta <- mod.mcml$x \%*\% beta + Z \%*\% Uk
        zeros <- rep(0, length(Uk))
        lfu[[k]] <- distRand(nu, Uk, mod.mcml$z, zeros)
        lfu.twid[[k]] <- distRand(nu.pql, Uk, mod.mcml$z, u.star)
        lfyu[[k]] <- el(mod.mcml$y, mod.mcml$x, eta, family.mcml)
        b[k] <- lfu[[k]]$value + lfyu[[k]]$value - lfu.twid[[k]]$value
    }
    a <- max(b)
    thing <- exp(b - a)
    value <- a - log(m) + log(sum(thing))
    v <- thing/sum(thing)
    Gpiece <- matrix(data = NA, nrow = nrow(umat), ncol = length(par))
    for (k in 1:nrow(umat)) {
        Gpiece[k, ] <- c(lfyu[[k]]$gradient, lfu[[k]]$gradient) * 
            v[k]
    }
    G <- apply(Gpiece, 2, sum)
    panda.list <- list()
    for (k in 1:nrow(umat)) {
        panda.list[[k]] <- c(lfyu[[k]]$gradient, lfu[[k]]$gradient) \%*\% 
            t(c(lfyu[[k]]$gradient, lfu[[k]]$gradient)) * v[[k]]
    }
    panda <- addMats(panda.list)
    lobster.list <- list()
    for (k in 1:nrow(umat)) {
        mat1 <- lfyu[[k]]$hessian
        mat2 <- lfu[[k]]$hessian
        d1 <- nrow(mat1)
        d2 <- nrow(mat2)
        newmat <- matrix(data = 0, nrow = d1 + d2, ncol = d1 + 
            d2)
        newmat[1:d1, 1:d1] <- mat1
        here <- d1 + 1
        there <- d1 + d2
        newmat[here:there, here:there] <- mat2
        lobster.list[[k]] <- newmat * v[k]
    }
    lobster <- addMats(lobster.list)
    hessian <- lobster + panda - G \%*\% t(G)
    list(value = value, gradient = G, hessian = hessian)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
