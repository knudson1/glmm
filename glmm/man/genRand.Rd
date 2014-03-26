\name{genRand}
\alias{genRand}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
genRand(sigma.star, s.star, z.list, m)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sigma.star}{
%%     ~~Describe \code{sigma.star} here~~
}
  \item{s.star}{
%%     ~~Describe \code{s.star} here~~
}
  \item{z.list}{
%%     ~~Describe \code{z.list} here~~
}
  \item{m}{
%%     ~~Describe \code{m} here~~
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
function (sigma.star, s.star, z.list, m) 
{
    nrand <- lapply(z.list, ncol)
    nrandom <- unlist(nrand)
    q <- sum(nrandom)
    if (q != length(s.star)) 
        stop("Number of random effects in dispute between s.star and mod.mcml$z")
    eek <- getEk(z.list)
    Aks <- Map("*", eek, sigma.star)
    A.star <- addVecs(Aks)
    u.star <- A.star * s.star
    u <- matrix(data = NA, nrow = m, ncol = q)
    for (k in 1:m) {
        u.swoop <- rnorm(q)
        u[k, ] <- u.swoop * A.star + u.star
    }
    list(u = u, u.star = u.star)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
