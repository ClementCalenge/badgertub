\name{nimbleToCoda}
\alias{nimbleToCoda}
\title{
  Converts the Result of a Nimble Fit to a Suitable Class for Use with
  Package coda
}
\description{
  This is used for the convergence diagnostics. Only the MCMC objects
  corresponding to our model can be used in this function.
}
\usage{
nimbleToCoda(sam, start = 1, end = nrow(sam[[1]]), thin = 1)
}
\arguments{
  \item{sam}{
    The MCMC samples returned by the function \code{nimbleMCMC} of the
    package nimble
  }
  \item{start}{
    the iteration number of the first observation
  }
  \item{end}{
    the iteration number of the last observation
  }
  \item{thin}{
    the thinning interval between consecutive observations
  }
}
\value{
  An object of class \code{mcmc.list}
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

## conversion for coda:
samdo <- (do.mcmc.out$samples)
ml <- nimbleToCoda(samdo, start=3001, end=1003000, thin=1000)


}
\keyword{model}
