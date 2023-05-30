\name{dpt}
\alias{dpt}
\docType{data}
\title{
  Map of the French Departments
}
\description{
  The object \code{dpt} is an object of class \code{"sf"} containing a
  map of the French departments (French administrative division).
}
\usage{data("dpt")}
\format{
  An object of class \code{"sf"} (see the package sf)
}
\source{
  Institut national de l'information geographique et forestiere. GEOFLA database.
}
\examples{
data(dpt)
plot(sf::st_geometry(dpt))
}
\keyword{datasets}
