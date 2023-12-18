simulateIndicator <-
    function(trappingPressure, cfc, intercept, rangeSlope=c(-0.4,0.4),
             tauu=0.73, rho = 0.04, nonadditive=FALSE, biasedSample=FALSE, nsim=1000,
             verbose=TRUE)
{

    ## The data.frame with commune/years
    expli <- expand.grid(ID=1:nrow(cfc), year=1:7)
    expli <- as.data.frame(expli)

    ## Sensitivity
    sens <- c(0.5, 0.5, 0.5, 0.75, 0.75, 0.75, 0.75)

    ## Preparation spatial matrix
    Hm <- spdep::listw2mat(spdep::nb2listw(spdep::poly2nb(as(cfc,"Spatial")), style="B"))
    diag(Hm) <- 1
    ni <- colSums(Hm)-1
    Q <- Hm
    Q[Q!=0] <- -1
    diag(Q) <- ni
    ei <- eigen(Q)
    ei$val <- 1/ei$val
    ei$vec <- ei$vec[,-ncol(ei$vec)]
    ei$val <- ei$val[-length(ei$val)]
    V <- ei$vec%*%diag(ei$val)%*%t(ei$vec)
    V <- V/tauu

    ## Distance for APC
    dis <- cfc  |> st_as_sfc() |> st_centroid()  |> st_coordinates()  |> dist()  |> as.matrix()
    dwij <- 1/(1+dis)

    ## Simulations
    lisi <- list()
    for (i in 1:nsim) {
        if (verbose)
            cat(i,"/",nsim,"\r")

        ## Simulation  spatial effects
        si <- MASS::mvrnorm(1,rep(0,ncol(V)),V)
        if (nonadditive) {
            si1 <- MASS::mvrnorm(1, rep(0, ncol(V)), V)
            si2 <- MASS::mvrnorm(1, rep(0, ncol(V)), V)
        }

        ## random SlopeYear
        slopeYear <- runif(1,rangeSlope[1],rangeSlope[2])

        ## Trapped animals
        if (biasedSample) {
            if (nonadditive) {
                bsi <- si1
            } else {
                bsi <- si
            }
            utm <- exp(bsi)
            utm <- utm/sum(utm)
            traPP <- trappingPressure*utm*length(bsi)
        } else {
            traPP <- trappingPressure
        }
        ssi <- MASS::rnegbin(nrow(expli), traPP, 0.48) ## using parameters estimated by glm.nb on the distri
        expli$Nb <- ssi

        ## model of prevalence
        if (nonadditive) {
            spatef <- sapply(1:nrow(expli), function(u) {
                (1-(expli$year[u]-1)/6)*si1[expli$ID[u]] + ((expli$year[u]-1)/6)*si2[expli$ID[u]]
            })
            si <- (si1+si2)/2
        } else {
            spatef <- si[expli$ID]
        }
        lpim <- intercept + spatef + expli$year * slopeYear
        pim <- exp(lpim)/(1+exp(lpim))
        lpimo <- intercept+si+4*slopeYear
        pimo <- exp(lpimo)/(1+exp(lpimo))
        pimo <- pimo[si>0]

        ## Simulation beta-binomiale
        a = pim * (1-rho)/rho
        b = (pim * rho - pim - rho + 1)/rho
        b2 <- (1-pim)*((1/rho)-1)
        p = rbeta(nrow(expli), a, b)
        y = rbinom(nrow(expli), expli$Nb, p*sens[expli$year])

        ## Calculation APC
        expli$Ninfect <- y
        expli$pr <- p
        expli$pim <- pim
        cpmsi <- apc4sim(expli, si, dwij)

        ## Indicator
        exb <- expli[expli$Nb>0,]
        exb$N0 <- exb$Nb-exb$Ninfect
        tab1 <- exb[rep(1:nrow(exb), exb$Ninfect),c("year","pr","ID")]
        tab0 <- exb[rep(1:nrow(exb), exb$N0),c("year","pr","ID")]
        tab1$y <- 1/sens[tab1$year]
        tab0$y <- 0
        tabt <- rbind(tab1,tab0)
        ## filter on high-prevalence commune
        tabt <- tabt[tabt$ID%in%c(1:length(si))[si>0],]
        explio <- expli[expli$ID%in%c(1:length(si))[si>0],]

        ## model
        tabt$year <- tabt$year-median(unique(tabt$year))
        mou <- lm(y~year, data=tabt)
        indic <- coef(mou)
        vc <- sqrt(diag(vcov(mou)))
        nt <- nrow(tabt)

        lisi[[i]] <- c(mean(pimo), cpmsi, indic, vc, nt, sum(si>0))
    }
    dfsi <- as.data.frame(do.call(rbind, lisi))
    names(dfsi) <- c("meanPrev","APC","InterceptReg","SlopeReg",
                     "SE_InterceptReg","SE_SlopeReg","Nind","Ninfect")

    return(dfsi)
}
