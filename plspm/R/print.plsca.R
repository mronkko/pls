`print.plsca` <-
function(x,...)
{
    Name <- c("$x.scores", "$x.wgs", "$x.loads", "$y.scores", "$y.wgs", "$y.loads",               
              "$cor.xt", "$cor.yu", "$cor.tu", "$cor.xu", "$cor.yt", "$R2X", "$R2Y", 
              "$com.xu", "$com.yt")
    Description <- c("X-scores (T)", "X-weights", "X-loadings", "Y-scores (U)", "Y-weights", 
         "Y-loadings", "X,T correlations", "Y,U correlations", "T,U correlations",
         "X,U correlations", "Y,T correlations", "explained variance of X by T",
         "explained variance of Y by U", "communality of X with U", "communality of Y with T")
    res1 <- cbind(Name, Description)
    cat("PARTIAL LEAST SQUARES CANONICAL ANALYSIS (PLS-CA)", "\n\n")
    cat("Results available in the following objects:", "\n")
    rownames(res1) <- 1:nrow(res1)
    print(res1, print.gap=3, justify="right")
    cat("---------------------------------------------------", "\n")    
    invisible(x)
}

