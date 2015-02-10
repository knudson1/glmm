\name{salamander}
\docType{data}
\alias{salamander}
\title{Salamander mating data set from McCullagh and Nelder (1989)}
\description{This data set presents the outcome of an experiment  conducted at the University of Chicago in 1986 to study interbreeding between populations of mountain dusky salamanders. }
 

\usage{data(salamander)}
\format{
  A data frame with the following columns:

   \describe{
\item{Mate}{Whether the salamanders mated (1) or did not mate (0).}
\item{Cross}{Cross between female and male type. A factor with four levels: \code{RR},\code{RW},\code{WR}, and \code{WW}. The type of the female salamander is listed first and the male is listed second. Rough Butt is represented by R and White Side is represented by W. For example, \code{Cross=WR} indicates a White Side female was crossed with a Rough Butt male.}
\item{Male}{Identification number of the male salamander. A factor.}
\item{Female}{Identification number of the female salamander. A factor.}

}
}


\references{McCullagh P. and Nelder, J. A. 1989. Generalized Linear Models, Section 14.5, Chapman and Hall/CRC.}

\examples{data(salamander)}
\keyword{generalized linear mixed models}

