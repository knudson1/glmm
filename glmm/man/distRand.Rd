\name{distRand}
\alias{distRand}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
distRand(nu, U, z.list, mu)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nu}{
%%     ~~Describe \code{nu} here~~
}
  \item{U}{
%%     ~~Describe \code{U} here~~
}
  \item{z.list}{
%%     ~~Describe \code{z.list} here~~
}
  \item{mu}{
%%     ~~Describe \code{mu} here~~
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
function (nu, U, z.list, mu) 
{
    T <- length(z.list)
    nrand <- lapply(z.list, ncol)
    nrandom <- unlist(nrand)
    totnrandom <- sum(nrandom)
    mu.list <- U.list <- NULL
    if (T == 1) {
        U.list[[1]] <- U
        mu.list[[1]] <- mu
    }
    if (T > 1) {
        U.list[[1]] <- U[1:nrandom[1]]
        mu.list[[1]] <- mu[1:nrandom[1]]
        for (t in 2:T) {
            thing1 <- sum(nrandom[1:t - 1]) + 1
            thing2 <- sum(nrandom[1:t])
            U.list[[t]] <- U[thing1:thing2]
            mu.list[[t]] <- mu[thing1:thing2]
        }
    }
    val <- gradient <- Hessian <- rep(0, T)
    for (t in 1:T) {
        you <- as.vector(U.list[[t]])
        mew <- as.vector(mu.list[[t]])
        Umu <- (you - mew) \%*\% (you - mew)
        val[t] <- as.numeric(-0.5 * nrandom[t] * log(nu[t]) - 
            Umu/(2 * nu[t]))
        gradient[t] <- -nrandom[t]/(2 * nu[t]) + Umu/(2 * (nu[t])^2)
        Hessian[t] <- nrandom[t]/(2 * (nu[t])^2) - Umu/((nu[t])^3)
    }
    value <- sum(val)
    if (T > 1) 
        hessian <- diag(Hessian)
    if (T == 1) 
        hessian <- matrix(Hessian, nrow = 1, ncol = 1)
    list(value = value, gradient = gradient, hessian = hessian)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
