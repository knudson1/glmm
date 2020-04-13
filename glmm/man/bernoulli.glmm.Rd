\name{bernoulli.glmm}
\alias{bernoulli.glmm}
\concept{generalized linear mixed model}
\concept{Monte Carlo}
\concept{likelihood approximation}

\title{Functions for the Bernoulli family. }
\description{Given a scalar \code{eta}, this calculates the cumulant and two derivatives for the Bernoulli family. Also checks that the data are entered correctly.}
\usage{
bernoulli.glmm()
}


\value{
\item{family.glmm }{The family name, as a string.}
\item{link}{The link function (canonical link is required), as a string.}
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
eta<--3:3
bernoulli.glmm()$family.glmm
bernoulli.glmm()$cum(eta)
bernoulli.glmm()$cp(1)
bernoulli.glmm()$cpp(2)

}
\keyword{models}
