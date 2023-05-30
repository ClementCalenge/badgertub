\name{giveXYhomoSens}
\alias{giveXYhomoSens}
\title{
  Converts a Sylvatub Dataset to a Dataframe and Attribute a Pair of
  Coordinates for each Animal.
}
\description{
  The function \code{giveXYhomoSens} converts a Sylvatub dataset (for a
  bTB cluster) into a dataframe with one row per animal, and give a pair
  of spatial coordinatess to each animal. These coordinates are randomly
  sampled in the commune where the animal was trapped. Then, the data
  are subsampled to allow a homogeneous sensitivity (see vignette of the
  package, appendix 4).
}
\usage{
giveXYhomoSens(tubData, tubConsts, cfc)
}
\arguments{
  \item{tubData}{
    The list of data used for the model fit in \code{nimbleMCMC}.
  }
  \item{tubConsts}{
    The list of constants used for the model fit in \code{nimbleMCMC}.
  }
  \item{cfc}{
    An object of class \code{sf} containing the map of the commune of
    the cluster.
  }
}
\value{
  A data.frame.
}
\author{
  Clement Calenge \email{clement.calenge@ofb.gouv.fr}
}
\seealso{
  \code{\link{nimbleData2df}}, which performs the same operation, but
  without attributing random coordinates to the animals.
}
\examples{

data(dotub)
data(cfcdo)

xydo <- giveXYhomoSens(dotub$data, dotub$consts, cfcdo)

plot(xydo[,c("X","Y")], asp=1)

}
\keyword{manip}
