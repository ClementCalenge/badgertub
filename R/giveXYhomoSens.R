giveXYhomoSens <-
function(tubData, tubConsts, cfc)
{
    dfxy <- nimbleData2df(tubData, tubConsts, keepHigh=FALSE)
    ta <- table(dfxy$ID)
    cfcr <- cfc[as.numeric(names(ta)),]
    sam <- st_sample(cfcr, ta)
    coor <- st_coordinates(sam)
    jo <- st_within(sam, cfc, sparse=FALSE)
    ID <- apply(jo,1,which)
    coor <- as.data.frame(coor)
    coor$ID <- ID
    dfxy <- dfxy[order(dfxy$ID),]
    coor <- coor[order(coor$ID),]
    if (!(all(ID==dfxy$ID)))
        stop("Problem with the function")
    dfxy$X <- coor[,1]
    dfxy$Y <- coor[,2]

    ## min sensitivity
    sensp <- min(dfxy$sensitivity)/dfxy$sensitivity
    dfxy$y <- rbinom(nrow(dfxy),size=dfxy$y, prob=sensp)
    return(dfxy)
}
