\name{radish2}
\docType{data}
\alias{radish2}
\concept{generalized linear mixed model}
\title{Radish count data set}
\description{Data on life history traits for the invasive California wild
     radish \emph{Raphanus sativus}. }

\usage{data(radish2)}
\format{
  A data frame with the following columns:

   \describe{
\item{Site}{Categorical. Experimental site where plant was grown. Two sites in this dataset.}
\item{Block}{Categorical. Blocked nested within site.}
\item{Region}{Categorical. Region from which individuals were obtained: northern, coastal California (N) or southern, inland California (S).}
\item{Pop}{Categorical. Wild population nested within region.}
\item{varb}{Categorical. Gives node of graphical model corresponding to each component of resp. This is useful for life history analysis (see aster package).}
\item{resp}{Response vector.}
\item{id}{Categorical. Indicates individual plants.}
}
}


\references{
% update to match update to help page in R package aster
These data are a subset of data previously analyzed using fixed effect
aster methods (R function \code{aster}) in the following.

Ridley, C. E. and Ellstrand, N. C. (2010).
Rapid evolution of morphology and adaptive life history in
    the invasive California wild radish (\emph{Raphanus sativus}) and
    the implications for management.
\emph{Evolutionary Applications}, \bold{3}, 64--76.

These data are a subset of data previously analyzed using random effect
aster methods (R function \code{reaster}) in the following.

Geyer, C. J., Ridley, C. E., Latta, R. G., Etterson, J. R.,
    and Shaw, R. G. (2013)
Local Adaptation and Genetic Effects on Fitness: Calculations
    for Exponential Family Models with Random Effects.
\emph{Annals of Applied Statistics}, \bold{7}, 1778--1795.
\doi{10.1214/13-AOAS653}.
}

\examples{data(radish2)}
\keyword{datasets}

