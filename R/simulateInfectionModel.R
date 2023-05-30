simulateInfectionModel <-
function(sam, tubConsts, tubData, verbose=TRUE)
{
    if (!is.list(sam))
        stop("sam should be a list")
    di <- unique(sapply(apply(sapply(sam, dim),1,unique, simplify=FALSE),length))
    if (any(di!=1))
        stop("sam should contain matrices of similar dimension")

    mcd <- do.call(rbind,sam)
    simy <- list()

    for (i in 1:nrow(mcd)) {
        if (verbose)
            cat(i, "/", nrow(mcd), "\r")
        intercept <- mcd[i,"intercept"]
        pente <- mcd[i,"slopeYear"]
        phi <- mcd[i,"phi"]
        si <- mcd[i,grep("^si", colnames(mcd))]
        tau <- mcd[i,ncol(mcd)]

        lpim <- intercept + si[tubConsts$ID] + tubData$year*pente
        pim <- exp(lpim)/(1+exp(lpim))

        a <- pim * phi
        b <- phi *(1-pim)

        pr <- rbeta(length(a), a, b)
        simy[[i]] <- rbinom(length(pr), prob=pr*tubData$sensitivity,size=tubData$Nb)
    }
    simy <- do.call(cbind,simy)
    return(simy)
}
