averageTrend <-
function(sam, tubData, tubConsts)
{

    samb <- do.call(rbind, sam)
    rne <- colMeans(samb)
    sime <- rne[grep("^si", names(rne))]
    quels <- sime>mean(sime)

    se <- apply(samb[,grep("^si",colnames(samb))], 2,
                function(x) x+samb[,"intercept"])
    sei <- se[,quels]
    y1 <- tubData$year
    y2 <- as.numeric(factor(tubData$year))
    moy <- y2[1]-y1[1]

    efan <- outer(samb[,"slopeYear"], sort(unique(y1)))
    rme <- sapply(1:nrow(efan), function(r) {
        g <- 1/(1+exp(-outer(sei[r,], efan[r,], "+")))
        mean(colMeans(apply(g,1,diff)))
    })

    rme <- rowMeans(apply(efan,2,function(x) {
        g <- apply(sei,2,function(y) 1/(1+exp(-x-y)))
    }))
    class(rme) <- c("APC","AL")
    return(rme)
}
