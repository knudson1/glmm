\name{bacteria}
\docType{data}
\alias{bacteria}
\title{Presence of Bacteria after Drug Treatments}
\description{Tests of the presence of the bacteria H. influenzae in children with otitis media in the Northern Territory of Australia.}

\usage{data(bacteria)}
\format{
  A data frame with the following columns:

   \describe{
    \item{y}{Presence or absence: a factor with levels \code{n} and \code{y}.}
    \item{ap}{active/placebo: a factor with levels \code{a} and \code{p}.}
    \item{hilo}{hi/low compliance: a factor with levels \code{hi} and \code{lo}.}
\item{week}{Numeric: week of test.}
\item{ID}{Subject ID: a factor.}
\item{trt}{A factor with levels \code{placebo}, \code{drug}, \code{drug+}, a re-coding of \code{ap} and \code{hilo}.}
\item{y2}{y reformatted as 0/1 rather than n/y.}
  }
}

\details{Dr. A. Leach tested the effects of a drug on 50 children with a history of otitis media in the Northern Territory of Australia. The children were randomized to the drug or the a placebo, and also to receive active encouragement to comply with taking the drug.

The presence of H. influenzae was checked at weeks 0, 2, 4, 6 and 11: 30 of the checks were missing and are not included in this data frame.}

\references{Venables, W. N. and Ripley, B. D. (2002) \emph{Modern Applied Statistics with S}, Fourth edition.  Springer.}
\examples{
data(bacteria)
}
\keyword{generalized linear mixed models}

