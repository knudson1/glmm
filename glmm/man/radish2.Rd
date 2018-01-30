\name{radish2}
\docType{data}
\alias{radish2}
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


\references{     These data are a subset of data previously analyzed using aster
     methods in the following.

     Ridley, C. E. and Ellstrand, N. C. (2010)  Rapid evolution of
     morphology and adaptive life history in the invasive California
     wild radish (\emph{Raphanus sativus}) and the implications for
     management.  \emph{Evolutionary Applications}, \bold{3}, 64--76.
}

\examples{data(radish2)}
\keyword{generalized linear mixed models}

