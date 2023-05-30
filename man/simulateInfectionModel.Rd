\name{simulateInfectionModel}
\alias{simulateInfectionModel}
\alias{compareSimActual}
\alias{print.csiac}
\title{
  Goodness-of-fit Checks on a Model
}
\description{
  The function \code{simulateInfectionModel} can be used to simulate
  fake dataset from a fitted model, to allow a check of the goodness of
  fit.  The function \code{compareSimActual} compares statistics
  calculated on the observed dataset with the credible intervals derived
  from the distribution of the same statistics calculated using
  simulated datasets.
}
\usage{
simulateInfectionModel(sam, tubConsts, tubData, verbose = TRUE)

compareSimActual(sim, tubData, tubConsts, levelCI=0.9)

\method{print}{csiac}(x, \dots)

}
\arguments{
  \item{sam}{
    The MCMC samples returned by the function \code{nimbleMCMC} of the
    package nimble.
  }
  \item{tubConsts}{
    The list of constants used for the model fit in \code{nimbleMCMC}.
  }
  \item{tubData}{
    The list of data used for the model fit in \code{nimbleMCMC}.
  }
  \item{verbose}{
    logical value indicating whether the function should return
    information on the progress of the simulation.
  }
  \item{sim}{
    Matrix returned by \code{simulateInfectionModel}.
  }
  \item{levelCI}{
    The level of the credible intervals used for the tests.
  }
  \item{x}{
    An object of class \code{"csiac"} returned by \code{compareSimActual}.
  }
  \item{\dots}{
    Additionnal elements to be passed to the generic function
    \code{print}.
  }
}
\value{
  \code{simulateInfectionModel} simulates a matrix with as many columns
  as there are samples in \code{sam}, and as many rows as there are
  commune-years in the datasets (in elements \code{y,Nb,sensitivity,}
  etc. of \code{tubData}.\cr

  \code{compareSimActual} returns a list of class \code{"csiac"} with
  elements named \code{total, peryear, percommune, percommuneyear,
  levelCI} containing the results of the comparison for the statistics
  "total number of infected animals",
  "number of infected animals for each year",
  "number of infected animals for each commune",
  "number of infected animals for each commune-year".  Elements named
  \code{total, peryear, percommune} are lists containing the
  observation(s) \code{obs}, the simulated values \code{sim}, the
  quantile \code{quant} of the distribution of the simulated values
  corresponding to the level \code{levelCI}, and the proportion of
  observed values within the limits of the credible intervals. The
  element \code{percommuneyear} contains only the elements \code{quant}
  and \code{prop}.
  
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


## simulation of a fake dataset for each sample (WARNING, VERY LONG --
## no need to execute, the dataset simDor loaded after these
## instructions already contains the results)
\dontrun{
  samdo <- (do.mcmc.out$samples)
  simDor <- simulateInfectionModel(samdo, dotub$consts, dotub$data)
}
data(simDor)


## Comparison here:
(csa <- compareSimActual(simDor, dotub$data, dotub$consts))






}
\keyword{model}
