summaryModel <-
function(sam, levelCI=0.9)
{
    if (!is.list(sam))
        stop("sam should be a list")
    di <- unique(sapply(apply(sapply(sam, dim),1,unique, simplify=FALSE),length))
    if (any(di!=1))
        stop("sam should contain matrices of similar dimension")

    samb <- do.call(rbind,sam)
    samb[,"phi"] <- 1/(1+samb[,"phi"])
    samb[,"tau"] <- 1/sqrt(samb[,"tau"])

    qs <- apply(samb,2,quantile,c((1-levelCI)/2,1-(1-levelCI)/2))
    ms <- apply(samb,2,median)
    ide <- c("slopeYear","phi", "tau")
    para <- ms[ide]
    parq <- apply(qs[,ide],2,
                  function(x)
        paste(" [", paste(round(x,2),collapse=", "), "]", sep=""))
    ob <- paste(round(para,2), parq)
    data.frame(Parameter=c("Slope of the year","Rho", "1/sqrt(tau)"),
               Estimates=ob)
}
