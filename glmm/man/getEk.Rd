\name{getEk}
\alias{getEk}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getEk(z)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{z}{
%%     ~~Describe \code{z} here~~
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
function (z) 
{
    nrand <- lapply(z, ncol)
    nrandom <- unlist(nrand)
    totnrandom <- sum(nrandom)
    nsigma <- length(z)
    Ekmat <- matrix(data = c(0), nrow = nsigma, ncol = totnrandom)
    Ek <- list()
    cumnrand <- c()
    cumnrand[1] <- 0
    for (i in 1:length(nrandom)) {
        cumnrand[i + 1] <- sum(nrandom[1:i])
    }
    for (i in 1:nsigma) {
        for (j in 1:totnrandom) {
            if (cumnrand[i] < j & cumnrand[i + 1] > j - 1) 
                Ekmat[i, j] <- 1
        }
        Ek[[i]] <- Ekmat[i, ]
    }
    Ek
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
