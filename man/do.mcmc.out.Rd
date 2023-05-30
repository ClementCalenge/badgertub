\name{do.mcmc.out}
\alias{do.mcmc.out}
\alias{bo.mcmc.out}
\alias{be.mcmc.out}
\docType{data}
\title{
  MCMC Samples From the Infection Model Fitted to the Sylvatub Dataset
  in the Three Main French Clusters
}
\description{
  The datasets \code{do.mcmc.out}, \code{bo.mcmc.out},
  \code{be.mcmc.out} contain the result of the fit by the package Nimble
  of the Infection Model of Calenge et al. (in prep.) to the data
  collected in the Dordogne/Charentes cluster, Burgundy cluster and
  Bearn cluster respectively (see examples to see how these datasets
  were obtained). 
}
\usage{

## Dordogne/Charentes
data("do.mcmc.out")

## Burgundy
data("bo.mcmc.out")

## Bearn
data("be.mcmc.out")

}
\format{
  These objects are returned by the function \code{nimbleMCMC} of the
  package \code{nimble}.
}
\references{

  Calenge C., Desvaux S., Payne A., Reveillaud E. and Richomme C. (in
  prep.) Assessing the bovine tuberculosis dynamics in three French badger
  populations with a participatory monitoring programme.

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


data(do.mcmc.out)
str(do.mcmc.out)


}
\keyword{datasets}
