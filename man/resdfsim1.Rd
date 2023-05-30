\name{resdfsim1}
\alias{resdfsim1}
\alias{resdfsim2}
\docType{data}
\title{
  Simulations Carried Out in Calenge et al. (in prep.)
}
\description{
  These two datasets contain the results of the simulations carried out
  in Calenge et al. (in prep.) to assess the suitability of the two
  indicators (mean prevalence level and mean prevalence trend) derived
  from the linear regression carried out on the Sylvatub dataset (see
  details and examples for how these datasets were obtained).
  
}
\usage{
data("resdfsim1")
}
\format{
  A data frame returned by the function \code{simulateIndicator} (see
  the help page of this function).
}
\details{
 The dataset \code{resdfsim1} corresponds to the first set of
 simulations carried out in Calenge et al. (in prep.) to estimate the
 trend over time of the infection (see examples).

 The dataset \code{resdfsim2} corresponds to the second set of
 simulations carried out in Calenge et al. (in prep.) to estimate the
 mean prevalence level in highly infected communes.
}
\references{

  Calenge C., Desvaux S., Payne A., Reveillaud E. and Richomme C. (in
  prep.) Assessing the bovine tuberculosis dynamics in three French badger
  populations with a participatory monitoring programme.
  
}
\examples{
#############################################
##
## First set of simulations

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


#############################################
##
## second set of simulations

## load the map of communes
data(cfcdo)

## Different intercepts (= different prevalence level)
inter <- -c(4:0)

## Different trapping pressures
trap <- c(0.5, 1, 3, 10)

## Simulations (WARNING, VERY SLOW, CAN TAKE SEVERAL HOURS)
\dontrun{
k <- 1
for (s in 1:4) {
  for (i in 1:5) {
    cat("** Prevalence level:", inter[i],"\n\n")
    sin <- simulateIndicator(trap[s], cfcdo, inter[i],
    rangeSlope=c(-0.4,0.4),
    nsim=1000, verbose=TRUE)
    sin$TrueIntercept <- inter[i]
    sin$TrapPress <- trap[s]
    res[[k]] <- sin
    k <- k+1
  }
}
resdfsim2 <- do.call(rbind, res)
}

data(resdfsim2)

## coverage
cocim <- coverageCI(resdfsim2$meanPrev,
                    resdfsim2$InterceptReg,
                    resdfsim2$SE_InterceptReg,
                    resdfsim2$Nind,
                    list(resdfsim2$TrapPress,
                         resdfsim2$TrueIntercept))
cocim$TrapPress <- paste0("mu = ", cocim$TrapPress)
names(cocim)[1] <- "Intercept"

cocim

}
\keyword{datasets}
