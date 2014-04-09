\name{addVecs}
\alias{addVecs}
\title{
Adds vectors of common length.
}
\description{Given a list of vectors of common length, this function sums component-wise and returns a single vector.}
\usage{addVecs(vecs)}

\arguments{
  \item{vecs}{A list of vectors that have equal length.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns a vector whose length is equal to the lengths of the incoming vectors.}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Christina Knudson}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
vec1<-1:6
veclist<-list(vec1,2*vec1,vec1+1)
addVecs(veclist)

## The function is currently defined as
function (vecs) 
{
    if (!is.list(vecs)) 
        vecs <- list(vecs)
    dvecs <- length(vecs[[1]])
    out <- rep(0, dvecs)
    for (d in 1:dvecs) {
        thing <- lapply(vecs, "[[", d)
        thing <- unlist(thing)
        out[d] <- sum(thing)
    }
    out
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ generalized linear mixed model }
\keyword{ Monte Carlo likelihood approximation }% __ONLY ONE__ keyword per line
