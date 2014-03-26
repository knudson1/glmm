\name{print.summary.mcla}
\alias{print.summary.mcla}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
print.summary.mcla(summ)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{summ}{
%%     ~~Describe \code{summ} here~~
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
function (summ) 
{
    stopifnot(inherits(summ, "summary.mcla"))
    cat("\nCall:\n", paste(deparse(summ$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")
    cat("Fixed Effects:")
    cat("\n")
    printCoefmat(summ$coefmat, digits = 3)
    cat("\n")
    cat("\n")
    cat("Variance Components for Random Effects (P-values are one-tailed):")
    cat("\n")
    printCoefmat(summ$nucoefmat, digits = 3)
    cat("\n")
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
