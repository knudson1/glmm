\name{cbpp2}
\docType{data}
\alias{cbpp2}
\title{Contagious Bovine Pleuropneumonia}
\description{This data set is a reformatted version of \code{cbpp} from the \code{lme4} package. Contagious bovine pleuropneumonia (CBPP) is a major disease of cattle in Africa, caused by a mycoplasma. This dataset describes the serological incidence of CBPP in zebu cattle during a follow-up survey implemented in 15 commercial herds located in the Boji district of Ethiopia. The goal of the survey was to study the within-herd spread of CBPP in newly infected herds. Blood samples were quarterly collected from all animals of these herds to determine their CBPP status. These data were used to compute the serological incidence of CBPP (new cases occurring during a given time period). Some data are missing (lost to follow-up). }

\usage{data(cbpp2)}
\format{
  A data frame with 3 columns:

   \describe{
    \item{Y}{Response vector. 1 if CBPP is observed, 0 otherwise.}
    \item{period}{A factor with levels 1 to 4.}
    \item{herd}{A factor identifying the herd (1 through 15).}

  }
}

\details{Serological status was determined using a competitive enzyme-linked immuno-sorbent assay (cELISA).}

\references{Lesnoff, M., Laval, G., Bonnet, P., Abdicho, S., Workalemahu, A., Kifle, D., Peyraud, A., Lancelot, R., Thiaucourt, F. (2004)  Within-herd spread of contagious bovine pleuropneumonia in Ethiopian highlands. \emph{Preventive Veterinary Medicine}, \bold{64}, 27--40. \doi{10.1016/j.prevetmed.2004.03.005}.}
\examples{
data(cbpp2)
}
\keyword{datasets}

