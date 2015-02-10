\name{summary.glmm}
\alias{summary.glmm}
\alias{print.summary.glmm}
\title{Summarizing GLMM Fits}
\description{ "summary" method for class \code{glmm} objects.}

\usage{
\method{summary}{glmm}(object, \dots)

\method{print}{summary.glmm}(x, digits = max(3, getOption("digits") - 3),
      signif.stars = getOption("show.signif.stars"), \dots)
}
\arguments{
  \item{object}{an object of class \code{glmm}, usually, resulting from a
    call to \code{\link{glmm}}.}
  \item{x}{an object of class \code{summary.glmm}, usually, a result of a
    call to \code{summary.glmm}.}
 \item{digits}{the number of significant digits to use when printing.}
  \item{signif.stars}{logical. If \code{TRUE}, \dQuote{significance stars}
    are printed for each coefficient.}
  \item{\dots}{further arguments passed to or from other methods.}

}

\value{
The function \code{summary.glmm} computes and returns a list of summary statistics of the fitted generalized linear mixed model given in \code{object}, using the components (list elements) \code{"call"} and \code{"terms"} from its argument, plus

 \item{coefficients}{a matrix for the fixed effects. The matrix has columns for the estimated coefficient, its standard error, t-statistic and corresponding (two-sided) p-value.}
 \item{nucoefmat}{a matrix with columns for the variance components. The matrix has columns for the estimated variance component, its standard error, t-statistic and corresponding (one-sided) p-value.}
\item{x}{the design matrix for the fixed effects.}
\item{z}{the design matrix for the random effects.}
\item{y}{the response vector.}
\item{fixedcall}{the call for the fixed effects.}
\item{randcall}{the call for the random effects.}
\item{family.mcml}{the family used to fit the model.}
\item{call}{the call to \code{\link{glmm}}.}
}

\author{Christina Knudson}




\seealso{The model fitting function \code{\link{glmm}}, the generic \code{\link{summary}}, and the function \code{\link{coef}}that extracts the fixed effect coefficients.
}


\keyword{generalized linear mixed models}
\keyword{models}
\keyword{Monte Carlo}
\keyword{maximum likelihood}
