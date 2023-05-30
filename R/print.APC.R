print.APC <-
function(x,...)
{
    if (!inherits(x, "APC"))
        stop("x should inherit class \"APC\"")
    re <- summary(x)
    if (inherits(x,"deltau")) {
        cat("*** Average predictive comparison for the year\n")
    } else {
        cat("*** Mean prevalence\n")
    }

    print(re)
}
