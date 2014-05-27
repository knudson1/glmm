\name{murder}
\docType{data}
\alias{murder}
\title{Number of Homicide Victims Known}
\description{Subjects responded to the question 'Within the past 12 months, how many people have you known personally that were victims of homicide?'}

\usage{data(murder)}
\format{
  A data frame with the following columns:

   \describe{
    \item{y}{The number of homicide victims known personally by the subject.}
    \item{race}{a factor with levels \code{black} and \code{white}.}
    \item{black}{a dummy variable to indicate whether the subject was black.}
\item{white}{a dummy variable to indiciate whether the subject was white.}
}
}

\references{Agresti, A. (2002) Categorical Data Analysis. Second edition. Wiley.}
\examples{
data(murder)
}
\keyword{generalized linear mixed models}

