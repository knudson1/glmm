\name{bernoulli.mcml}
\alias{bernoulli.mcml}

\title{
Functions for the Bernoulli family. }
\description{
Given a scalar eta, this calculates the cumulant and two derivatives for the Bernoulli family. Also checks that the data are entered correctly.}
\usage{
bernoulli.mcml()
}

\details{
 }
\value{
\item{family.mcml }{The family name, as a string.}
  \item{cum }{The cumulant function.}
  \item{cp }{The first derivative of the cumulant function.}
  \item{cpp }{The second derivative of the cumulant function.}
  \item{checkData }{A function to check that all data are either 0 or 1.}
}
\references{}
\author{Christina Knudson}
\note{This function is to be used by the mcml command.}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{mcml}}
}
\examples{
bernoulli.mcml()$family.mcml
bernoulli.mcml()$cum(2)
bernoulli.mcml()$cp(2)
bernoulli.mcml()$cpp(2)

}

\keyword{ generalized linear mixed model }
\keyword{ Monte Carlo likelihood approximation }
