`print.plsreg2` <-
function(x, ...)
{
    Name <- c("$x.scores", "$x.loads", "$y.scores", "$y.loads", "$raw.wgs", "$weights", 
              "$cor.tx", "$cor.ty", "$std.coef", "$coeffs", "$y.pred", "$resid",
              "$expvar", "$Q2", "$Q2cum", "$VIP")
    Description <- c("X-scores", "X-loadings", "Y-scores", "Y-loadings", "raw weights",
           "modified weights", "correlations between X-T", "correlations between Y-T", 
           "standard coefficients", "coefficients", "Y-predicted", "residuals",
           "explained variance", "Q2 index", "cummulated Q2", "variable importance for projection")
    res1 <- cbind(Name, Description)
    cat("PARTIAL LEAST SQUARES REGRESSION 2 (PLS-R2)", "\n\n")
    cat("Results available in the following objects:", "\n")
    rownames(res1) <- 1:nrow(res1)
    print(res1, print.gap=3, justify="right")
    cat("---------------------------------------------------", "\n")    
    invisible(x)
}

