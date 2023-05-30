compareSimActual <-
function(sim, tubData, tubConsts, levelCI=0.9)
{

    ## Total
    tots <- colSums(sim)
    toto <- sum(tubData$y)
    qtot <- quantile(tots, c((1-levelCI)/2, 1-(1-levelCI)/2))
    qtot <- paste0("[",qtot[1],",",qtot[2],"]")

    ## For each year
    ans <- sapply(sort(unique(tubData$year)),function(i) colSums(sim[tubData$year==i,]))
    qpa <- apply(ans,2,quantile,c((1-levelCI)/2, 1-(1-levelCI)/2))
    ano <- tapply(tubData$y, tubData$year, sum)
    paqo <- mean(ano>=qpa[1,]&ano<=qpa[2,])

    ## For each commune
    parans <- sapply(sort(unique(tubConsts$ID)), function(i) colSums(sim[tubConsts$ID==i,,drop=FALSE]))
    qa <- apply(parans,2,quantile,c((1-levelCI)/2, 1-(1-levelCI)/2))
    parano <- tapply(tubData$y, tubConsts$ID, sum)
    pqo <- mean(parano>=qa[1,]&parano<=qa[2,])

    ## For each year and commune
    ict <- apply(sim,1,quantile,c((1-levelCI)/2, 1-(1-levelCI)/2))
    qi <- mean(tubData$y>=ict[1,]&tubData$y<=ict[2,])

    ## return results
    lieval <- list(total=list(obs=toto, sim=tots, quant=qtot, prop=1),
                   peryear=list(obs=ano, sim=ans, quant=qpa, prop=paqo),
                   percommune=list(obs=parano, sim=parans, quant=qa, prop=pqo),
                   percommuneyear=list(quaant=ict, prop=qi),
                   levelCI=levelCI)
    class(lieval) <- "csiac"
    return(lieval)
}
