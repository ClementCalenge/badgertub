nimbleToCoda <-
function(sam, start=1, end=nrow(sam[[1]]), thin=1)
{
    if (!is.list(sam))
        stop("sam should be a list")
    di <- unique(sapply(apply(sapply(sam, dim),1,unique, simplify=FALSE),length))
    if (any(di!=1))
        stop("sam should contain matrices of similar dimension")

    ml <- mcmc.list(lapply(sam, function(chu)
    {
        r <- cbind(chu[,c("phi","slopeYear","tau")],
                   do.call(cbind,
                           lapply(grep("^si", colnames(chu)),
                                  function(i) chu[,"intercept"]+chu[,i])))
        colnames(r) <- c(c("phi","slopeYear","tau"),
                         colnames(chu)[grep("^si", colnames(chu))])
        mcmc(r, start=start, end=end, thin=thin)
    }))
    return(ml)
}
