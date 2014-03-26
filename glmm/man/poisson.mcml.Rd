\name{poisson.mcml}
\alias{poisson.mcml}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
poisson.mcml()
}
%- maybe also 'usage' for other objects documented here.
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
function () 
{
    family.mcml <- "poisson.mcml"
    cum <- function(eta) exp(eta)
    cp <- function(eta) exp(eta)
    cpp <- function(eta) exp(eta)
    checkData <- function(x) {
        bads <- sum(x < 0)
        if (bads > 0) 
            stop("data must be nonnegative.")
        return(NULL)
    }
    out <- list(family.mcml = family.mcml, cum = cum, cp = cp, 
        cpp = cpp, checkData = checkData)
    class(out) <- "mcml.family"
    return(out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
