apc4sim <-
function(expli, si, dwij)
{

    explib <- expli[expli$Nb>0,]
    explib <- explib[explib$ID%in%c(1:length(si))[si>0],]

    dpred <- as.vector(outer(explib$pr,explib$pr,"-"))
    du <- as.vector(outer(explib$year, explib$year, "-"))
    nn <- as.vector(outer(explib$Nb, explib$Nb))
    wij <- unlist(lapply(1:nrow(explib), function(i) dwij[explib$ID,
                                                          explib$ID[i]]))

    resu <- sum(nn*wij*dpred*sign(du))/sum(nn*wij*du*sign(du))
    return(resu)
}
