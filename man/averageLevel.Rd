\name{averageLevel}
\alias{averageLevel}
\alias{averageTrend}
\title{
  Calculate the Average Prevalence Level and the Average Trend in
  Prevalence Using the Simple Regression Approach.
}
\description{
  These two functions calculate the average prevalence level and trends
  using the simple regression approach advocated in Calenge et al. (in
  prep.). 
}
\usage{

averageLevel(sam, tubData, tubConsts)

averageTrend(sam, tubData, tubConsts)

}
\arguments{
  \item{sam}{
    The MCMC samples returned by the function \code{nimbleMCMC} of the
    package nimble.
  }
  \item{tubData}{
    The list of data used for the model fit in \code{nimbleMCMC}.
  }
  \item{tubConsts}{
    The list of constants used for the model fit in \code{nimbleMCMC}.   
  }
}
\value{
  A vector of class \code{"APC"} with one level/trend estimated per MCMC sample.
}
\author{
  Clement Calenge \email{clement.calenge@ofb.gouv.fr}
}
\seealso{
  \code{\link{apc}}
}
\examples{

## Load the dataset
data(dotub)

## Load the results of the model fit 
data(do.mcmc.out)
samdo <- (do.mcmc.out$samples)

## Level and trend
(aldo <- averageLevel(samdo,dotub$data, dotub$consts))
(atdo <- averageTrend(samdo,dotub$data, dotub$consts))

}
\keyword{model}
