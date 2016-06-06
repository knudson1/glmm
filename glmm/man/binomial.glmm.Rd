\name{binomial.glmm}
\alias{binomial.glmm}

\title{Functions for the Binomial family. }
\description{Given a scalar \code{eta} and the number of trials, this calculates the cumulant and two derivatives for the Bernoulli family. Also checks that the data are entered correctly.}
\usage{
bernoulli.glmm()
}


\value{
\item{family.glmm }{The family name, as a string.}
\item{link}{The link function (canonical link is required), as a string.}
  \item{cum }{The cumulant function.}
  \item{cp }{The first derivative of the cumulant function.}
  \item{cpp }{The second derivative of the cumulant function.}
  \item{checkData }{A function to check that all data are nonnegative.}
}

\author{Christina Knudson}
\note{This function is to be used by the \code{\link{glmm}} command.}



\seealso{
\code{\link{glmm}}
}
\examples{
eta<--3:3
ntrials <- 1
binomial.glmm()$family.glmm
binomial.glmm()$cum(eta, ntrials)
binomial.glmm()$cp(1, ntrials)
binomial.glmm()$cpp(2, ntrials)

}

\keyword{ generalized linear mixed model }
\keyword{ Monte Carlo likelihood approximation }
