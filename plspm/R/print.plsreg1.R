`print.plsreg1` <-
function(x, ...)
{
    Name <- c("$scores", "$x.loads", "$y.loads", "$u.scores", "$raw.wgs", "$mod.wgs", 
              "$std.coef", "$coeffs", "$R2", "$y.pred", "$resid", "$cor.sco", "$T2")
    Description <- c("X-scores", "X-loadings", "Y-loadings", "U-scores", "raw weights",
         "modified weights", "standard coefficients", "coefficients", "R-squared", 
         "y-predicted", "residuals", "score correlations", "T2 hotelling")
    if (length(x)==15) {
        Name <- c(Name, "$Q2")
        Description <- c(Description, "Q2 cross validation")
    }
    res1 <- cbind(Name, Description)
    cat("PARTIAL LEAST SQUARES REGRESSION 1 (PLS-R1)", "\n\n")
    cat("Results available in the following objects:", "\n")
    rownames(res1) <- 1:nrow(res1)
    print(res1, print.gap=3, justify="right")
    cat("---------------------------------------------------", "\n")    
    invisible(x)
}

