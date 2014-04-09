\name{catVecs}
\alias{catVecs}

\title{
Vector concatenation
}
\description{
Concatenates vectors into one vector.
}
\usage{
catVecs(vecs)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{vecs}{A list of vectors of any length.}
}
\details{
The vectors are concatenated in the order given.
}
\value{Returns the vector that results from concatenating the list of vectors vecs.}
\references{
%% ~put references to the literature/web site here ~
}
\author{Christina Knudson}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
vec1<-c(1:6)
vec2<-c(10:18)
vec3<-c(4,7)
catVecs(list(vec1,vec2,vec3))

## The function is currently defined as
function (vecs) 
{
    if (!is.list(vecs)) 
        vecs <- list(vecs)
    lengths <- lapply(vecs, length)
    lengths <- unlist(lengths)
    out <- rep(0, sum(lengths))
    starthere <- 1
    for (d in 1:length(vecs)) {
        endhere <- starthere + lengths[d] - 1
        out[c(starthere:endhere)] <- unlist(vecs[[d]])
        starthere <- starthere + lengths[d]
    }
    out
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{generalized linear mixed model}
\keyword{ monte carlo likelihood approximation } 
