print.csiac <-
function(x, ...)
{
    if (!inherits(x,"csiac"))
        stop("x should be of class csiac")
    lev <- x$levelCI
    cat("**********\n** Goodness of fit of the model\n")
    cat("\nTotal number of infected animals:\n","obs =", x$total$obs,
        "   Expected",paste0(lev*100, "% CI ="), x$total$quant,"\n")
    cat("\nNumber of infected animals per year\n",
        "Proportion of years in the", paste0(lev*100, "% CI:"), x$peryear$prop,"\n")
    cat("\nNumber of infected animals per commune\n",
        "Proportion of commune in the", paste0(lev*100, "% CI:"), x$percommune$prop,"\n")
    cat("\nNumber of infected animals per year-commune\n",
        "Proportion of year-commune in the", paste0(lev*100, "% CI:"), x$percommuneyear$prop,"\n")
}
