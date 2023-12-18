## ----setup, include=FALSE, cache=FALSE----------------------------------------
# set global chunk options
library('knitr')
opts_chunk$set(fig.path="badgertub-",
               fig.align="center",
               fig.show="hold",
               echo=TRUE,
               results="markup",
               fig.width=10,
               fig.height=10, out.width='\\linewidth',
               out.height='\\linewidth',
               cache=FALSE,
               dev='png',
               concordance=TRUE,
               error=FALSE)
opts_knit$set(aliases = c(h = 'fig.height',
              w = 'fig.width',
              wo='out.width',
              ho='out.height'))
options(replace.assign=TRUE,width=60)
set.seed(9567)


## ----eval=FALSE---------------------------------------------------------------
## ## If devtools is not yet installed, type
## install.packages("devtools")
## 
## ## Install the package badgertub
## devtools::install_github("ClementCalenge/badgertub", ref="main")


## ----load-badgertub-----------------------------------------------------------
library(badgertub)


## ----load-datasets------------------------------------------------------------
## Two maps
data(dpt)
data(cfcdo)

## The Sylvatub dataset
data(dotub)

str(dotub)


## ----load-nimble--------------------------------------------------------------
library(nimble)


## ----code-nimble-infection-model----------------------------------------------
Tubcode <- nimbleCode({
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


## ----starting-values-model----------------------------------------------------
dotubInits <- list(tau = 0.1, slopeYear = 0.2,
                   intercept=-3, phi=20,
                   si=rnorm(dotub$consts$P),
                   p=rep(0.1, dotub$consts$N))


## ----model-fit, eval=FALSE----------------------------------------------------
## set.seed(777)
## do.mcmc.out <- nimbleMCMC(code = Tubcode, constants = dotub$consts,
##                           data = dotub$data, inits = dotubInits,
##                           nchains = 4, niter = 1003000, nburnin=3000, thin=1000,
##                           summary = TRUE, WAIC = TRUE,
##                           monitors = c('tau',"intercept","slopeYear","phi", "si"))


## ----load-results-model-fit---------------------------------------------------
data(do.mcmc.out)


## ----conversion-results-to-coda-----------------------------------------------
samdo <- (do.mcmc.out$samples)
ml <- nimbleToCoda(samdo, start=3001, end=1003000, thin=1000)


## ----plot-chains, eval=FALSE--------------------------------------------------
## plot(ml)


## ----gelman-diag-complete, eval=FALSE-----------------------------------------
## library(coda)
## gelman.diag(ml, multivariate=FALSE)


## ----load-coda-silent, echo=FALSE---------------------------------------------
library(coda)


## ----gelman-diag-top-parameters-----------------------------------------------
mlb <- ml
for (i in 1:4)
    mlb[[i]] <- mlb[[i]][,1:3]
gelman.diag(mlb,multivariate=FALSE)


## ----simulate-goodness-of-fit, eval=FALSE-------------------------------------
## simDor <- simulateInfectionModel(samdo, dotub$consts, dotub$data)


## ----load-simus-gof-----------------------------------------------------------
data(simDor)


## ----compare-sim-actual-gof---------------------------------------------------
(csa <- compareSimActual(simDor, dotub$data, dotub$consts))


## ----estimated-parameters-of-the-model----------------------------------------
summaryModel(samdo)


## ----map-spatial-random-effects-----------------------------------------------
showSpatialEffects(samdo,cfcdo,dpt)


## ----calculate-APC, eval=FALSE------------------------------------------------
## APCdo <- apc(samdo, cfcdo, dotub$consts, dotub$data)


## ----load-result-calculation-APC----------------------------------------------
data(APCdo)
APCdo


## ----true-mean-prevalence-level-----------------------------------------------
aldo <- averageLevel(samdo,dotub$data, dotub$consts)
aldo


## ----indicators-dordogne------------------------------------------------------
dfb <- nimbleData2df(dotub$data, dotub$consts, sam=samdo)
dfb$year <- dfb$year-4
dfb$response <- dfb$y/dfb$sensitivity
me <- lm(response~year, data=dfb)
summary(me)


## -----------------------------------------------------------------------------
## A list of length two (one element per simulated situation).  Each
## element is a vector containing the limits of the range within which
## the slope of the year is randomly sampled:
##
## - Situation 1: the slope of the year is increasing (i.e. positive
## slope), and is randomly sampled between 0 and 0.4)
##
## - Situation 2: the slope is either increasing or decreasing, and is
## - sampled between -0.4 and 0.4.
lirs <- list(c(0,0.4),
             c(-0.4,0.4))
##
## Four trapping pressures are defined here (mean number of trapped
## animals per commune)
trap <- c(0.5, 1, 3, 10)

## Two intercepts for the model (either low prevalence = -3.1 or high
## prevalence -1.38).
situ <- c(-3.1, -1.38)


## ----first-set-simulations, eval=FALSE----------------------------------------
## k <- 1
## res <- list()
## for (i in 1:4) {
##     cat("** Trapping pressure:", trap[i],"\n\n")
##     for (s in 1:2) {
##         cat("Situation:", s,"\n\n")
##         sin <- simulateIndicator(trap[i], cfcdo, situ[s],
##                                  rangeSlope=lirs[[s]],
##                                  nsim=1000, verbose=TRUE)
##         sin$Situation <- c("Low Increasing","High")[s]
##         sin$TrapPress <- paste0("mu = ",trap[i])
##         res[[k]] <- sin
##         k <- k+1
##     }
## }
## resdfsim1 <- do.call(rbind, res)
## resdfsim1$TrapPress <- factor(resdfsim1$TrapPress,
##                               levels=c("mu = 0.5", "mu = 1",
##                                        "mu = 3", "mu = 10"))


## -----------------------------------------------------------------------------
data(resdfsim1)
head(resdfsim1)


## ----plot-first-set-simus-fig2-paper------------------------------------------
library(ggplot2)
ggplot(resdfsim1, aes(x=APC, y=SlopeReg))+
    geom_point(aes(col=TrapPress), alpha=0.5)+
    geom_abline(slope=1,intercept=0, lwd=1)+
    facet_grid(Situation~TrapPress)+coord_fixed()+
    xlab("True prop. of animals becoming infected in one year")+
    ylab("Estimated value by the regression")+
    theme(legend.position = "none")


## ----estimation-mean-prevalence-first-set-simulations-------------------------
ggplot(resdfsim1, aes(x=meanPrev, y=InterceptReg))+
    geom_abline(slope=1,intercept=0, lwd=1, col="red")+
    geom_point(alpha=0.5)+
    facet_grid(Situation~TrapPress)+coord_fixed()+
    xlab("True mean Prevalence")+
    ylab("Estimated value by the regression")


## ----coverage-CI-first-set-of-simulations-------------------------------------
cocia <- coverageCI(resdfsim1$APC,
                   resdfsim1$SlopeReg,
                   resdfsim1$SE_SlopeReg,
                   resdfsim1$Nind,
                   list(resdfsim1$TrapPress,
                        resdfsim1$Situation))
cocia


## ----parameters-second-set-simulations----------------------------------------
## The intercepts
inter <- -c(4:0)

## The simulated trapping pressures
trap <- c(0.5, 1, 3, 10)


## ----second-set-of-simulations, eval=FALSE------------------------------------
## k <- 1
## for (s in 1:4) {
##     for (i in 1:5) {
##         cat("** Prevalence level:", inter[i],"\n\n")
##         sin <- simulateIndicator(trap[s], cfcdo, inter[i],
##                                  rangeSlope=c(-0.4,0.4),
##                                  nsim=1000, verbose=TRUE)
##         sin$TrueIntercept <- inter[i]
##         sin$TrapPress <- trap[s]
##         res[[k]] <- sin
##         k <- k+1
##     }
## }
## resdfsim2 <- do.call(rbind, res)
## resdfsim2$TrapPress <- paste0("mu = ", resdfsim2$TrapPress)
## resdfsim2$TrapPress <- factor(resdfsim2$TrapPress,
##                               levels=c("mu = 0.5", "mu = 1",
##                                        "mu = 3", "mu = 10"))
## resdfsim2$TrueIntercept <- paste0("alpha = ", resdfsim2$TrueIntercept)
## resdfsim2$TrueIntercept <- factor(resdfsim2$TrueIntercept,
##                                   levels=c("alpha = -4", "alpha = -3",
##                                            "alpha = -2", "alpha = -1",
##                                            "alpha = 0"))


## ----load-results-second-set-of-simulations-----------------------------------
data(resdfsim2)
head(resdfsim2)


## ----plot-results-second-set-of-simulations-----------------------------------
ggplot(resdfsim2, aes(x=meanPrev, y=InterceptReg))+
    geom_abline(slope=1,intercept=0, lwd=1, col="red")+
    geom_point(alpha=0.5)+
    facet_grid(TrueIntercept~TrapPress)+coord_fixed()+
    xlab("True mean Prevalence")+
    ylab("Estimated value by the regression")


## ----FOI-second-set-simus-APC-------------------------------------------------
ggplot(resdfsim2, aes(x=APC, y=SlopeReg))+
    geom_abline(slope=1,intercept=0, lwd=1, col="red")+
    geom_point(alpha=0.5)+
    facet_grid(TrueIntercept~TrapPress)+coord_fixed()+
    xlab("True prop. of animals becoming infected in one year")+
    ylab("Estimated value by the regression")


## ----coverage-proba-of-CI-on-mean-prev----------------------------------------
cocim <- coverageCI(resdfsim2$meanPrev,
                    resdfsim2$InterceptReg,
                    resdfsim2$SE_InterceptReg,
                    resdfsim2$Nind,
                    list(resdfsim2$TrapPress,
                         resdfsim2$TrueIntercept))
names(cocim)[1] <- "Intercept"
cocim


## ----third-set-of-simulations, eval=FALSE-------------------------------------
## res <- list()
## interb <- c(0,-2)
## k <- 1
## for (b in 1:6) {
##     for (i in 1:2) {
##         cat("** Intercept:", i,"\n\n")
## 
##         sin <- simulateIndicator(2, cfcdo, interb[i],
##                                  rangeSlope=c(-0.4,0.4),
##                                  biasedSample=((b%%2)==0),
##                                  nonadditive=(b>2&b<5),
##                                  tau=c(0.1,0.73)[(b>2)+1],
##                                  nsim=1000, verbose=TRUE)
##         sin$TrueIntercept <- interb[i]
##         sin$biasedSample <- ((b%%2)==0)
##         sin$nonadditive <- (b>2&b<5)
##         sin$tau <- c(0.1,0.73)[(b>2)+1]
##         res[[k]] <- sin
##         k <- k+1
##     }
## }
## 
## resdfsim3 <- do.call(rbind, res)
## resdfsim3$Situation <- paste0(c("Additive","Interaction")[resdfsim3$nonadditive+1],", ",
##                               c("moderate","strong")[(resdfsim3$tau<0.2)+1])
## resdfsim3$Sampling <- c("random","directed")[resdfsim3$biasedSample+1]
## resdfsim3$TrueIntercept <- paste0("alpha = ", resdfsim3$TrueIntercept)
## resdfsim3$TrueIntercept <- factor(resdfsim3$TrueIntercept,
##                                   levels=c("alpha = -2",
##                                            "alpha = 0"))
## 


## -----------------------------------------------------------------------------
data(resdfsim3)
head(resdfsim3)


## ----FOI-third-set-simus-APC--------------------------------------------------
ggplot(resdfsim3, aes(x=APC, y=SlopeReg))+
    geom_abline(slope=1,intercept=0, lwd=1, col="red")+
    geom_point(alpha=0.5)+
    facet_grid(Sampling~Situation)+coord_fixed()+
    xlab("APC")+
    ylab("Estimated slope by the regression")



## ----mean-third-set-simus-APC-------------------------------------------------
ggplot(resdfsim3, aes(x=meanPrev, y=InterceptReg))+
    geom_abline(slope=1,intercept=0, lwd=1, col="red")+
    geom_point(alpha=0.5)+
    facet_grid(Sampling~Situation)+coord_fixed()+
    xlab("Mean prevalence")+
    ylab("Estimated intercept by the regression")



## ----coverage-proba-of-CI-on-mean-prev-3--------------------------------------
cocim3 <- coverageCI(resdfsim3$meanPrev,
                    resdfsim3$InterceptReg,
                    resdfsim3$SE_InterceptReg,
                    resdfsim3$Nind,
                    list(resdfsim3$Situation,
                         resdfsim3$Sampling))
names(cocim3)[1] <- "Intercept"
cocim3


## ----coverage-proba-of-CI-on-slope-3------------------------------------------
cocis3 <- coverageCI(resdfsim3$APC,
                    resdfsim3$SlopeReg,
                    resdfsim3$SE_SlopeReg,
                    resdfsim3$Nind,
                    list(resdfsim3$Situation,
                         resdfsim3$Sampling))
names(cocis3)[1] <- "Slope"
cocis3


## ----give-random-coordinates-Dordogne-----------------------------------------
set.seed(777)
xydo <- giveXYhomoSens(dotub$data, dotub$consts, cfcdo)
head(xydo)


## ----GAM-Dordogne-Charentes---------------------------------------------------
library(mgcv)
## Full model with space-time interactions and smoother for the year
gamdo_full <- gam(y~s(X,Y, k=40)+s(year,k=3)+
                      ti(X,Y,year, d=c(2,1), k=c(30,3)),
                  data=xydo, family=binomial, method="REML")


## ----summary-full-model-------------------------------------------------------
summary(gamdo_full)


## ----give-random-coordinates-Burgundy-----------------------------------------
data(botub)
data(cfcbo)

set.seed(777)
xybo <- giveXYhomoSens(botub$data, botub$consts, cfcbo)


## ----GAM-Burgundy-------------------------------------------------------------
## Full model with space-time interactions and smoother for the year
gambo_full <- gam(y~s(X,Y, k=40)+s(year,k=3)+
                      ti(X,Y,year, d=c(2,1), k=c(30,3)),
                  data=xybo, family=binomial, method="REML")


## ----summary-full-model-Burgundy----------------------------------------------
summary(gambo_full)


## ----give-random-coordinates-Bearn--------------------------------------------
data(betub)
data(cfcbe)

set.seed(777)
xybe <- giveXYhomoSens(betub$data, betub$consts, cfcbe)


## ----GAM-Bearn----------------------------------------------------------------
## Full model with space-time interactions and smoother for the year
gambe_full <- gam(y~s(X,Y, k=40)+s(year,k=3)+
                      ti(X,Y,year, d=c(2,1), k=c(30,3)),
                  data=xybe, family=binomial, method="REML")


## ----summary-full-model-Bearn-------------------------------------------------
summary(gambe_full)

