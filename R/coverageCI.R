coverageCI <-
function(trueval, estimated, se, n, listf, levelCI=0.95)
{
    if (length(listf)!=2)
        stop("listf should be of length 2")
    if (length(unique(c(length(trueval), length(estimated),
                        length(se), length(n),
                        length(listf[[1]]), length(listf[[2]]))))>1)
        stop("non-convenient data")
    ci <- t(sapply(1:length(trueval),
                   function(i)
        estimated[i]+c(-1,1)*qt(1-(1-levelCI)/2, n[i]-2)*se[i]))
    coverint <- (tapply(ci[,1]<trueval&ci[,2]>=trueval,
                        listf, mean))
    class(coverint) <- "table"
    coverint <- as.data.frame(coverint)
    coverint[,1:2] <- coverint[,2:1]
    names(coverint) <- c("Situation","TrapPress", "Probability")
    return(coverint)
}
