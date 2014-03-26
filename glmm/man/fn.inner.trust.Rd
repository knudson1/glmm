\name{fn.inner.trust}
\alias{fn.inner.trust}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
fn.inner.trust(mypar, Y, X, Z, A, family.mcml, nbeta)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mypar}{
%%     ~~Describe \code{mypar} here~~
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
  \item{A}{
%%     ~~Describe \code{A} here~~
}
  \item{family.mcml}{
%%     ~~Describe \code{family.mcml} here~~
}
  \item{nbeta}{
%%     ~~Describe \code{nbeta} here~~
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
function (mypar, Y, X, Z, A, family.mcml, nbeta) 
{
    beta <- mypar[1:nbeta]
    s <- mypar[-(1:nbeta)]
    eta <- X \%*\% beta + Z \%*\% A \%*\% s
    family.mcml <- getFamily(family.mcml)
    mu <- family.mcml$cp(eta)
    value <- el(Y, X, eta, family.mcml)$value - 0.5 * s \%*\% s
    db <- t(X) \%*\% (Y - mu)
    ds <- A \%*\% t(Z) \%*\% (Y - mu) - s
    gradient <- c(db, ds)
    cdub <- family.mcml$cpp(eta)
    cdub <- as.vector(cdub)
    cdub <- diag(cdub)
    kyoo <- nrow(A)
    piece1 <- (-t(X) \%*\% cdub \%*\% X)
    piece2 <- (-t(X) \%*\% cdub \%*\% Z \%*\% A)
    piece3 <- (-A \%*\% t(Z) \%*\% cdub \%*\% X)
    piece4 <- (-A \%*\% t(Z) \%*\% cdub \%*\% Z \%*\% A - diag(kyoo))
    hessian <- rbind(cbind(piece1, piece2), cbind(piece3, piece4))
    list(value = value, gradient = gradient, hessian = hessian)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
