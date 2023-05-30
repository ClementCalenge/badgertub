\name{simulateIndicator}
\alias{simulateIndicator}
\alias{coverageCI}
\title{
  Simulate a Sylvatub Dataset and Calculate Prevalence Trends and Levels Indicators.
}
\description{
  The function \code{simulateIndicator} can be used to simulate the
  infection model by Calenge et al. (in prep.) and estimate the
  prevalence level and trend indicators derived from the simple
  regression level (see details). The function
  \code{coverageCI} can be used to estimate the coverage probability of
  the confidence intervals derived for the simulated indicators.
}
\usage{

simulateIndicator(trappingPressure, cfc, intercept,
                  rangeSlope = c(-0.4, 0.4), nsim = 1000,
                  verbose = TRUE)

coverageCI(trueval, estimated, se, n, listf, levelCI = 0.95)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{trappingPressure}{
    The mean number of animals to simulate per commune.
  }
  \item{cfc}{
    An object of class \code{sf} containing the map of the commune of
    the cluster.
  }
  \item{intercept}{
    The intercept in the model of prevalence.
  }
  \item{rangeSlope}{
    The range of slopes for the year.
  }
  \item{nsim}{
    The number of simulations.
  }
  \item{verbose}{
    logical value indicating whether the function should return
    information on the progress of the simulation.
  }
  \item{trueval}{
    A vector of true values for the parameter for which the coverage
    probability is to be assessed.
}
  \item{estimated}{
    A vector of estimated values for the parameter for which the coverage
    probability is to be assessed.
  }
  \item{se}{
    A vector of standard errors corresponding to the estimated values
    for the parameter for which the coverage probability is to be
    assessed. 
  }
  \item{n}{
    The number of individuals used in the linear regression used to
    derive the indicators.
  }
  \item{listf}{
    A list of two factors indicating the levels for which the coverage
    probability is to be calculated (e.g. factors indicating different
    trapping pressures and different epidemiological situations, see
    examples). 
  }
  \item{levelCI}{
    The level of the confidence intervals used for the tests.
  }
}
\details{
  
  This function (i) simulates a Sylvatub dataset in the
  Dordogne/Charentes cluster by drawing a number of tested animals from
  a beta-binomial distribution with mean \code{trappingPressure} and
  dispersion \code{theta=0.48}); (ii) simulates for each tested animal a
  result of the test using the Bayesian model developped in Calenge et
  al. (in prep.), using the conditional autoregressive precision
  \code{tau = 0.73} for the commune effects, \code{intercept} as
  intercept of the model (see Calenge et al. in prep.), and a slope for
  the year randomly sampled in the range defined by \code{rangeSlope},
  (iii) calculates the average predictive comparison, mean prevalence
  level and prevalence level and trends indicators derived from the
  simple linear regression described by Calenge et al. (in prep.).
  
}
\value{
  
  A data frame with one row per simulation and containing the estimated
  mean prevalence and APC in highly infected communes, the value of the
  true intercept and slope of the simple linear regression approach
  (indicators), the standard error of these estimated values, the number
  of simulated individuals and the number of simulated infected
  individuals.

}
\references{
  
  Calenge C., Desvaux S., Payne A., Reveillaud E. and Richomme C. (in
  prep.) Assessing the bovine tuberculosis dynamics in three French badger
  populations with a participatory monitoring programme.

}
\author{
    Clement Calenge \email{clement.calenge@ofb.gouv.fr}
}
\seealso{
  \code{\link{apc}}
}
\examples{

## Load the map of communes
data(cfcdo)

## Two sets of simulations for the slope of the year (two ranges) 
lirs <- list(c(0,0.4),c(-0.4,0.4))

## Four tested trapping pressure
trap <- c(0.5, 1, 3, 10)

## Two intercepts
situ <- c(-3.1, -1.38)

## WARNING!!! Very slow (can take several hours to execute). The
## result resdfsim1 is available as a dataframe in the package (loaded
## after this part).
\dontrun{
k <- 1
res <- list()
for (i in 1:4) {
    cat("** Trapping pressure:", trap[i],"\n\n")
    for (s in 1:2) {
        cat("Situation:", s,"\n\n")
        sin <- simulateIndicator(trap[i], cfcdo, situ[s],
                                 rangeSlope=lirs[[s]],
                                 nsim=1000, verbose=TRUE)
        sin$Situation <- c("Low Increasing","High")[s]
        sin$TrapPress <- paste0("mu = ",trap[i])
        res[[k]] <- sin
        k <- k+1
    }
}
resdfsim1 <- do.call(rbind, res)

resdfsim1$TrapPress <- factor(resdfsim1$TrapPress,
                              levels=c("mu = 0.5", "mu = 1", "mu = 3", "mu = 10")) 
}

## Load the resulting dataset
data(resdfsim1)

## plot the results
ggplot2::ggplot(resdfsim1, ggplot2::aes(x=APC, y=SlopeReg))+
    ggplot2::geom_point(ggplot2::aes(col=TrapPress), alpha=0.5)+
    ggplot2::geom_abline(slope=1,intercept=0, lwd=1)+
    ggplot2::facet_grid(Situation~TrapPress)+ggplot2::coord_fixed()+
    ggplot2::xlab("True prop. of animals becoming infected in one year")+
    ggplot2::ylab("Estimated value by the regression")+
    ggplot2::theme(legend.position = "none")
    
## The coverage probability:
cocia <- coverageCI(resdfsim1$APC,
                   resdfsim1$SlopeReg,
                   resdfsim1$SE_SlopeReg,
                   resdfsim1$Nind,
                   list(resdfsim1$TrapPress,
                   resdfsim1$Situation))
cocia


}
\keyword{model}
