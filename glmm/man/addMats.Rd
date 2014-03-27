\name{addMats}
\alias{addMats}
\title{Adds a list of matrices with common dimensions.}
\description{Performs element-wise addition for any number of matrices in a list.}
\usage{addMats(matList)
}
\arguments{
  \item{matList}{A list of matrices that all have identical dimensions.}
}
\details{The number of rows/columns must be identical for each matrix in matList.}
\value{Returns a matrix whose dimensions are equal to the dimensions of the incoming matrices.}
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
mat1<-matrix(data=1:6,nrow=2)
matlist<-list(mat1,2*mat1)
addMats(matlist)

## The function is currently defined as
function (matList) 
{
    if (!is.list(matList)) 
        matList <- list(matList)
    ncells <- length(as.vector(matList[[1]]))
    nr <- nrow(matList[[1]])
    matAsVec <- lapply(matList, as.vector)
    out <- rep(0, ncells)
    for (d in 1:ncells) {
        thing <- lapply(matAsVec, "[[", d)
        thing <- unlist(thing)
        out[d] <- sum(thing)
    }
    matrix(data = out, nrow = nr)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ generalized linear mixed model}
\keyword{ Monte Carlo likelihood approximation}
