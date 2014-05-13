\name{poisson.glmm}
\alias{poisson.glmm}

\title{Functions for the Poisson family.}
\description{Given a scalar \code{eta}, this calculates the cumulant and two derivatives for the Poisson family. Also checks that the data are entered correctly.}
\usage{poisson.glmm()}


\value{
\item{family.glmm }{The family name, as a string.}
  \item{cum }{The cumulant function.}
  \item{cp }{The first derivative of the cumulant function.}
  \item{cpp }{The second derivative of the cumulant function.}
  \item{checkData }{A function to check that all data are either 0 or 1.}
}

\author{Christina Knudson}

\note{This function is to be used by the \code{\link{glmm}} command.}




\seealso{
\code{\link{glmm}}
}
\examples{
poisson.glmm()$family.glmm
poisson.glmm()$cum(2)
poisson.glmm()$cp(2)
poisson.glmm()$cpp(2)
}

\keyword{ generalized linear mixed model }
\keyword{ Monte Carlo likelihood approximation }
