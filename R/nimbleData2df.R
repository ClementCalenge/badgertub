nimbleData2df <-
function(tubData,tubConsts, sam=NA, keepHigh=TRUE)
{
    if (!is.list(tubData))
        stop("non-convenient data")
    if (!is.list(tubConsts))
        stop("non-convenient data")
    if (!all(c("adj", "num", "weights", "Nb",
               "y", "sensitivity", "year")%in%names(tubData)))
        stop("non-convenient data")
    if (!all(c("N", "P", "L", "ID")%in%names(tubConsts)))
        stop("non-convenient data")
    if (keepHigh&!is.list(sam))
        stop("sam cannot be missing when keepHigh is TRUE")

    if (is.list(sam)) {
        di <- unique(sapply(apply(sapply(sam, dim),1,unique, simplify=FALSE),length))
        if (any(di!=1))
            stop("sam should contain matrices of similar dimension")
        samdo <- do.call(rbind,sam)
        com <- colMeans(samdo)
        sirr <- com["intercept"]+com[grep("^si", names(com))]
        quels <- sirr>mean(sirr)
    }
    re <- as.data.frame(tubData[c("y","sensitivity","year","Nb")])
    re$ID <- tubConsts$ID
    re2 <- do.call(rbind,
                   lapply(1:nrow(re),
                          function(i) {
                       res <- re[rep(i,re$Nb[i]),]
                       res$y <- c(rep(1,res$y[1]),
                                  rep(0,res$Nb[1]-res$y[1]))
                       return(res)
        }))
    if (keepHigh) {
        re2 <- re2[re2$ID%in%c(1:length(quels))[quels],]
    }
    re2$year <- as.numeric(factor(re2$year))
    return(re2)
}
