print.plspm.groups <-
function(x,...)
{
    cat("GROUP COMPARISON IN PLS-PM FOR PATH COEFFICIENTS", "\n\n")
    cat("Scale of Data:", "\t\t", x$settings[[1]], "\n")
    cat("Weighting Scheme:", "\t", x$settings[[2]], "\n")
    cat("Selected method: ", "\t", x$settings[[3]], "\n")
    if (x$settings[[3]]=="bootstrap") {
        cat("Number of resamples: ", "\t", x$reps, "\n\n")
    } else
        cat("Number of permutations: ", x$reps, "\n\n")
    cat("$test", "\n")
    print(x$test, print.gap=2)
    cat("\n")
    cat("Inner models in the following objects:", "\n")
    cat("$global ", "\n")
    cat("$group1 ", "\n")
    cat("$group2 ", "\n")
}

