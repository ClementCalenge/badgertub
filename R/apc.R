apc <-
function(sam, cfc, tubConsts, tubData, nsim=nrow(sam[[1]])*length(sam), verbose=TRUE)
{
    if (!is.list(sam))
        stop("sam should be a list")
    di <- unique(sapply(apply(sapply(sam, dim),1,unique, simplify=FALSE),length))
    if (any(di!=1))
        stop("sam should contain matrices of similar dimension")
    if (!inherits(cfc,"sf"))
        stop("cfc should inherit the class sf")

    samdo <- do.call(rbind,sam)
    samdo[,"phi"] <- 1/(1+samdo[,"phi"])
    samdo[,"tau"] <- 1/sqrt(samdo[,"tau"])

    com <- colMeans(samdo)

    lisi <- 0
    sirr <- com["intercept"]+com[grep("^si", names(com))]
    quels <- sirr>mean(sirr)

    ## Random sumpling of the iterations
    sam <- sample(1:nrow(samdo), nsim)

    for (r in 1:nsim) {
        if (verbose)
            cat(r,"/", nsim, "\r")
        i <- sam[r]

        ## recuup des effets spatiaux
        si <- samdo[i,grep("^si", names(com))]
        intercept <- samdo[i,"intercept"]
        slopeYear <- samdo[i,"slopeYear"]

        lpim <- intercept+si[tubConsts$ID]+tubData$year*slopeYear
        pim <- exp(lpim)/(1+exp(lpim))

        ## Simulation beta-binomiale
        rho = samdo[i,"phi"]
        a = pim * (1-rho)/rho
        b = (pim * rho - pim - rho + 1)/rho
        b2 <- (1-pim)*((1/rho)-1)
        p = rbeta(length(tubData$y), a, b)

        dfa <- data.frame(ID=tubConsts$ID,
                          pr=p,
                          year=tubData$year,
                          Nb=tubData$Nb)


        dis <- cfc  |> st_as_sfc() |> st_centroid()  |> st_coordinates()  |> dist()  |> as.matrix()
        dwij <- 1/(1+dis)

        ## Application de la formule:
        cpmsi <- capado(dfa, quels, dwij)
        lisi[r] <- cpmsi
    }
    class(lisi) <- c("APC","deltau")
    return(lisi)
}
