showSpatialEffects <-
function(sam, cfc, dpt, rangeCol=c(0,0.3))
{
    if (!is.list(sam))
        stop("sam should be a list")
    di <- unique(sapply(apply(sapply(sam, dim),1,unique, simplify=FALSE),length))
    if (any(di!=1))
        stop("sam should contain matrices of similar dimension")
    if (!inherits(dpt,"sf"))
        stop("dpt should inherit the class sf")
    if (!inherits(cfc,"sf"))
        stop("cfc should inherit the class sf")

    samb <- do.call(rbind,sam)
    sir <- apply(samb[,grep("^si",colnames(samb))], 2, function(x) mean(x+samb[,"intercept"]))
    cfcos <- cfc
    cfcos$zeta <- exp(sir)/(1+exp(sir))
    bb <- st_bbox(cfc)
    araj <- max(c(diff(bb[c(1,3)]), diff(bb[c(2,4)])))/6
    dor<-ggplot() + geom_sf(data = cfcos, aes(fill = .data$zeta), col=NA)+
        geom_sf(data=dpt, fill=NA, size=1, col="black")+
        geom_sf(data=dpt, fill=NA, size=0.5, col="green")+
        scale_fill_viridis(option="plasma",limits=rangeCol,name="Prevalence")+
        coord_sf(xlim = c(bb[1]-araj,bb[3]+araj),
                 ylim = c(bb[2]-araj,bb[4]+araj), expand = FALSE)+
        theme_void()+theme(panel.border = element_rect(colour = "black", fill=NA, size=2))
    return(dor)
}
