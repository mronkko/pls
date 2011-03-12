print.plspm <-
function(x, ...)
{
    k <- length(x) - 2
    Name <- paste(rep("$",k),names(x)[1:k], sep="")
    Description <- c("outer model", "inner model", "scaled LVs", "LVs for scaled=FALSE",
                     "outer weights", "loadings", "path coefficients matrix", "R-squared")
    if (k==13) # no bootstrap
    {
        if (is.null(x$data)) {
            Description <- c(Description, "outer correlations", "summary inner model", "total effects",
                         "unidimensionality", "goodnes-of-fit")
        } else {
            Name <- c(Name, "$data")
            Description <- c(Description, "outer correlations", "summary inner model", "total effects",
                         "unidimensionality", "goodnes-of-fit", "data matrix")
        }
    }
    if (k==14) # with bootstrap
    {
        if (is.null(x$data)) {
    Description <- c(Description, "outer correlations", "summary inner model", "total effects",
                        "unidimensionality", "goodnes-of-fit", "bootstrap results")
        } else {
            Name <- c(Name, "$data")
    Description <- c(Description, "outer correlations", "summary inner model", "total effects",
                        "unidimensionality", "goodnes-of-fit", "bootstrap results", "data matrix")
        }
    }             
    res1 <- cbind(Name, Description)
    cat("\n")
    cat("PARTIAL LEAST SQUARES PATH MODELING (PLS-PM)", "\n")
    if (k==8) cat("             Basic Algorithm", "\n")               
    cat("----------------------------------------------", "\n")    
    cat("Results available in the following objects:", "\n\n")
    rownames(res1) <- 1:nrow(res1)
    print(res1, print.gap=3, justify="right")
    cat("----------------------------------------------", "\n")
    cat("You can also use the function 'summary'", "\n\n")    
    invisible(x)
}

