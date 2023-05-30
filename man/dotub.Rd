\name{dotub}
\alias{dotub}
\alias{botub}
\alias{betub}
\docType{data}
\title{
  Bovine Tuberculosis Dataset Collected by the Sylvatub Network
}
\description{
  The datasets \code{dotub}, \code{botub} and \code{betub} contain the
  data collected by the Sylvatub network on the badger in the
  Dordogne/Charentes cluster, Burgundy cluster and Bearn cluster
  respectively. 
}
\usage{data("dotub")}
\format{
  A List containing two named elements: the element named \code{data}
  contains the data used for the model fit and the elements named
  \code{consts} contains the constants used in this fit.  More precisely,
  the element \code{data} contains:
  \describe{
    \item{adj}{vector of indices of the adjacent locations (neighbors)
      of each spatial location.  This is a sparse representation of the
      full adjacency matrix. Used in \code{dcar_normal} in the package
      \code{nimble}.}
    \item{num}{vector giving the number of neighboring locations of each
      spatial location, with length equal to the total number of
      locations.  Used in \code{dcar_normal} in the package
      \code{nimble}.}
    \item{weights}{vector of symmetric unnormalized weights associated
      with each pair of adjacent locations, of the same length as
      \code{adj}.  Also used in \code{dcar_normal} in the package
      \code{nimble}.}
    \item{Nb}{For a given commune and a given year, the number of
      trapped badgers that were analyzed by the Sylvatub network.}
    \item{y}{For a given commune and a given year, the number of
      trapped badgers that were analyzed by the Sylvatub network and
      were identified as positive.}
    \item{sensitivity}{For a given commune and a given year, the
      sensitivity of the tests used (either 0.5 for a bacterial culture,
      or 0.75 for a PCR test).}
    \item{year}{The centred year used in the model.}
  }
  The element \code{consts} contains:
  \describe{
    \item{N}{The number of commune-years in \code{data} (number of
      elements of \code{y,Nb,sensitivity}, etc.)}
    \item{P}{The number of communes of the cluster}
    \item{L}{The number of elements in \code{adj}}
    \item{ID}{The commune ID for each element in the data
      (\code{y,Nb,}etc.)}
  }
}
\details{
  These are the original data.
}
\source{
  Sylvatub Network.
}
\examples{
## no examples
}
\keyword{datasets}
