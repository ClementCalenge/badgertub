capado <-
function(df, quels, dwij)
{

    explib <- df[df$ID%in%c(1:length(quels))[quels],]

    dpred <- as.vector(outer(explib$pr,explib$pr,"-"))
    du <- as.vector(outer(explib$year, explib$year, "-"))
    nn <- as.vector(outer(explib$Nb, explib$Nb))
    wij <- unlist(lapply(1:nrow(explib), function(i) dwij[explib$ID,
                                                          explib$ID[i]]))

    resu <- sum(nn*wij*dpred*sign(du))/sum(nn*wij*du*sign(du))
    return(resu)
}
