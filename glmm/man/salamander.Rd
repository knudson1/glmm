\name{salamander}
\docType{data}
\alias{salamander}
\title{Salamander mating data set from McCullagh and Nelder (1989)}
\description{This data set presents the outcome of an experiment  conducted at the University of Chicago in 1986 to study interbreeding between populations of mountain dusky salamanders. }
 

\usage{data(salamander)}
\format{
  A data frame with the following columns:

   \describe{
    \item{Season}{A factor with levels \code{Spring} and \code{Summer} to categorize when the experiment was conducted.}
    \item{Experiment}{A factor with levels \code{1},\code{2}, \code{3} to indicate the experiment number.}
    \item{TypeM}{Type of the male salamander; Rough Butt=R and White Side=W.}
    \item{TypeF}{Type of the female salamander; Rough Butt=R and White Side=W.}
\item{Cross}{Cross between female and male type. A factor with four levels: \code{RR},\code{RW},\code{WR}, and \code{WW}. For example, \code{Cross=WR} indicates a White Side female was crossed with a Rough Butt male.}
\item{Male}{Identification number of the male salamander. A factor.}
\item{Female}{Identification number of the female salamander. A factor.}
\item{Mate}{Whether the salamanders mated (1) or did not mate (0).}
}
}


\references{McCullagh P. and Nelder, J. A. 1989. Generalized Linear Models, Section 14.5, Chapman and Hall/CRC.}

\examples{data(salamander)}
\keyword{generalized linear mixed models}

