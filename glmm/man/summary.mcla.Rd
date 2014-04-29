\name{summary.mcla}
\alias{summary.mcla}
\alias{print.summary.mcla}
\title{Summarizing GLMM Fits}
\description{ "summary" method for class \code{mcla} objects.}
\usage{summary.mcla(object)}
\arguments{
  \item{object}{an object of class \code{"mcla"}, usually, resulting from a
    call to \code{\link{glmm}}.
}
}

\value{
The function \code{summary.mcla} computes and returns a list of summary statistics of the fitted generalized linear mixed model given in \code{object}, using the components (list elements) \code{"call"} and \code{"terms"} from its argument, plus

 \item{coefficients}{a matrix for the fixed effects. The matrix has columns for the estimated coefficient, its standard error, t-statistic and corresponding (two-sided) p-value.}
 \item{nucoefmat}{a matrix with columns for the variance components. The matrix has columns for the estimated variance component, its standard error, t-statistic and corresponding (one-sided) p-value.}

}

\author{Christina Knudson}




\seealso{The model fitting function\code{\link{glmm}}, \code{\link{summary}}.

Function \code{\link{coef}} will extract the fixed effect coefficients.
}


\keyword{generalized linear mixed models}
\keyword{models}
\keyword{Monte Carlo}
\keyword{maximum likelihood}
