\name{cfcdo}
\alias{cfcdo}
\alias{cfcbo}
\alias{cfcbe}
\docType{data}
\title{
  Maps of the Communes in the Three bovine Tuberculosis Clusters in France.
}
\description{
  The datasets \code{cfcdo, cfcbo, cfcbe} contain the maps (class
  \code{"sf"} in the package sf) of the communes of the bTB clusters.
}
\usage{
data("cfcdo")

data("cfcbo")

data("cfcbe")

}
\format{
  A data frame of class \code{"sf"}.
}
\source{
  IGN France.
}
\examples{

data(cfcdo)
plot(cfcdo)

## And similarly for other clusters

}
\keyword{datasets}
