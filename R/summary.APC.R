summary.APC <-
function(object, ...)
{
    if (!inherits(object, "APC"))
        stop("object should inherit class \"APC\"")
    res <- c(mean(object), sd(object))
    if (inherits(object,"deltau")) {
        names(res) <- c("Delta_u", "SE")
    } else {
        names(res) <- c("MeanPrev","SE")
    }
    return(res)
}
