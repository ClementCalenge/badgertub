\name{summaryModel}
\alias{summaryModel}
\alias{showSpatialEffects}
\title{
  Main Results from a Model Fitted to a Sylvatub Dataset
}
\description{
  \code{summaryModel} summarize the fit of a model to the Sylvatub
  Dataset. \code{showSpatialEffects} plots a map of the spatial
  effects.
}
\usage{
summaryModel(sam, levelCI = 0.9)

showSpatialEffects(sam, cfc, dpt, rangeCol=c(0,0.3))
}
\arguments{
  \item{sam}{
    The MCMC samples returned by the function \code{nimbleMCMC} of the
    package nimble.
  }
  \item{levelCI}{
    The levels of the credible interval in the summary.
  }
  \item{cfc}{
    An object of class \code{sf} containing a map of the communes of the
    cluster under study.
  }
  \item{dpt}{
    An object of class \code{sf} containing a map of the departments of
    the cluster under study.
  }
  \item{rangeCol}{
    A numeric vector of length two providing limits of
    the scale. Use NA to refer to the existing minimum
    or maximum (passed to \code{scale_fill_viridis}).
  }
}
\value{
  \code{summaryModel} returns a data.frame with the summary.
  \code{showSpatialEffects} returns a ggplot.
}
\author{
  Clement Calenge \email{clement.calenge@ofb.gouv.fr}
}
\examples{

## Load data of the Dordogne/Charentes cluster
data(dotub)

## Load the results of the model fit 
data(do.mcmc.out)

## (for the record, the fitted model was):
Tubcode <- nimble::nimbleCode({
    tau ~ dgamma(1,1)
    slopeYear ~ dnorm(0, 0.01)
    intercept ~ dnorm(0, 0.01)
    phi ~ dgamma(0.1,0.1)

    si[1:P] ~ dcar_normal(adj[1:L], weights[1:L], num[1:P], tau)

    for (i in 1:N) {
        lpim[i] <- intercept + si[ID[i]] + year[i]*slopeYear
        pim[i] <- exp(lpim[i])/(1+exp(lpim[i]))

        a[i] <- pim[i] * phi
        b[i] <- phi *(1-pim[i])

        p[i] ~ dbeta(a[i], b[i])
        y[i] ~ dbinom(prob=p[i]*sensitivity[i],size=Nb[i])
    }
})

## and is fitted by (WARNING, VERY LONG -- no need to execute, the
## dataset do.mcmc.out loaded before already contains the results)
\dontrun{
dotubInits <- list(tau = 0.1, slopeYear = 0.2,
                   intercept=-3, phi=20,
                   si=rnorm(dotub$consts$P),
                   p=rep(0.1, dotub$consts$N))

do.mcmc.out <- nimble::nimbleMCMC(code = Tubcode,
                                  constants = dotub$consts,
                                  data = dotub$data,
                                  inits = dotubInits, nchains = 4,
                                  niter = 1003000, nburnin=3000,
                                  thin=1000, summary = TRUE,
                                  WAIC = TRUE,
                                  monitors = c("tau","intercept","slopeYear","phi", "si")))
}


## Summary of the model
samdo <- (do.mcmc.out$samples)
summaryModel(samdo)

## Plot of the spatial effects
data(cfcdo)
data(dpt)
showSpatialEffects(samdo,cfcdo,dpt)

}
\keyword{model}
