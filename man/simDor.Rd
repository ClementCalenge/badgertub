\name{simDor}
\alias{simDor}
\alias{simBou}
\alias{simBea}
\docType{data}
\title{
  Check of the Goodness of Fit of the Infection Models Fitted to the
  Sylvatub Datasets in the Three French Clusters.
}
\description{
  The datasets \code{simDor}, \code{simBou} and \code{simBea} are
  simulated Sylvatub datasets obtained from the models to the data
  collected in the Charentes/Dordogne cluster, Burgundy cluster and
  Bearn cluster respectively, using the function
  \code{simulateInfectionModel} (see examples below for the code that
  was used to generate them). The function \code{compareSimActual} can
  be used on them to assess the goodness of fit of these models
}
\usage{

data("simDor")

data("simBou")

data("simBea")

}
\format{
  A matrix returned by the function \code{simulateInfectionModel}.
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


## Simulate the model (also very slow, resulting data in simDor dataset, loaded after)
samdo <- (do.mcmc.out$samples)
simDor <- simulateInfectionModel(samdo, dotub$consts, dotub$data)

}

## We obtain
data(simDor)
head(simDor)

## And we assess the goodness of fit (see ?compareSimActual)
(csa <- compareSimActual(simDor, dotub$data, dotub$consts))



}
\keyword{datasets}
