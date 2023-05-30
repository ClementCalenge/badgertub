\name{apc}
\alias{apc}
\alias{summary.APC}
\alias{print.APC}
\title{
  Average Predictive Comparison for the Different Models of the bovine
  Tuberculosis Infection in the Three French Clusters.
}
\description{
  The function \code{apc} calculates the average predictive comparison
  for the year in a model of the bovine Tuberculosis infection in a
  cluster. 
}
\usage{
apc(sam, cfc, tubConsts, tubData,
    nsim = nrow(sam[[1]]) * length(sam), verbose = TRUE)

\method{summary}{APC}(object, \dots)

\method{print}{APC}(x, \dots)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sam}{
    The MCMC samples returned by the function \code{nimbleMCMC} of the
    package nimble.
  }
  \item{cfc}{
    An object of class \code{sf} containing a map of the communes of the
    cluster under study.
  }
  \item{tubConsts}{
    The list of constants used for the model fit in \code{nimbleMCMC}.
  }
  \item{tubData}{
    The list of data used for the model fit in \code{nimbleMCMC}.
  }
  \item{nsim}{
    The number of simulations used to estimate the APC.
  }
  \item{verbose}{
    Logical value indicating whether the function should return
    information on the progress of the simulation.
  }
  \item{object,x}{
    Object of class \code{"APC"} returned by the function \code{apc}.
  }
  \item{\dots}{
    Additional arguments passed to the generic functions.
  }
}
\value{
  The function \code{apc} returns a list of class \code{"APC"}
  containing the APC calculated for each sample in \code{sam}.
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

do.mcmc.out <- nimbleMCMC(code = Tubcode, constants = dotub$consts,
                          data = dotub$data, inits = dotubInits,
                          nchains = 4, niter = 1003000, nburnin=3000, thin=1000,
                          summary = TRUE, WAIC = TRUE,
                          monitors = c("tau","intercept","slopeYear","phi", "si")))
}


## simulation of a fake dataset for each sample (WARNING, VERY LONG --
## no need to execute, the dataset APCdo loaded before already
## contains the results)
data(cfcdo)
\dontrun{
  samdo <- (do.mcmc.out$samples)
  APCdo <- apc(samdo, cfcdo, dotub$consts, dotub$data)
}
data(APCdo)
APCdo


}
\keyword{model}
