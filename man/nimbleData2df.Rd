\name{nimbleData2df}
\alias{nimbleData2df}
\title{
  Converts a Sylvatub Dataset to a Dataframe
}
\description{
  \code{nimbleData2df} converts a Sylvatub dataset formatted for Nimble
  into a dataframe with one row per animal and with the columns
  corresponding to the variables in the model (y = infected/not
  infected, sensitivity of the test, year, Nb of animals in the same
  commune during the same year).
}
\usage{
nimbleData2df(tubData, tubConsts, sam = NA, keepHigh = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tubData}{
    The list of data used for the model fit in \code{nimbleMCMC}.
  }
  \item{tubConsts}{
    The list of constants used for the model fit in \code{nimbleMCMC}.
  }
  \item{sam}{
    Optionnally, the MCMC samples returned by the function
    \code{nimbleMCMC} of the package nimble (only used when
    \code{keepHigh=TRUE}). 
}
  \item{keepHigh}{
    Logical value indicating whether only animals from highly infected
    communes (i.e. with a spatial effect greater than average) should be
    included in the resulting dataframe.
  }
}
\value{
  A dataframe.
}
\author{
  Clement Calenge \email{clement.calenge@ofb.gouv.fr}
}
\examples{

## Load data of the Dordogne/Charentes cluster
data(dotub)

## Load the results of the model fit 
data(do.mcmc.out)
samdo <- (do.mcmc.out$samples)

## converts to data frame
dfb <- nimbleData2df(dotub$data, dotub$consts, sam=samdo)



}
\keyword{manip}
