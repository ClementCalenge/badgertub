averageLevel <-
function(sam, tubData, tubConsts)
{

    samb <- do.call(rbind, sam)
    samb2 <- samb[,"intercept"]+samb[,grep("^si", colnames(samb))]
    sime <- colMeans(samb2)
    un <- sort(unique(tubConsts$ID))
    sime <- sime[un]

    ## Highly infected communes
    quels <- sime>mean(sime)


    ## middle year
    med <- median(unique(tubData$year))

    se <- samb[,grep("^si",colnames(samb))]+samb[,"intercept"]+samb[,"slopeYear"]*med
    se <- se[,un]
    sei <- se[,quels]

    rme <- rowMeans(1/(1+exp(-sei)))

    class(rme) <- c("APC","AL")
    return(rme)
}
