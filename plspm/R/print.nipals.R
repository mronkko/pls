`print.nipals` <-
function(x, ...)
{
    Name <- c("$values", "$scores", "$loadings", "$cor.sco", 
              "$disto", "$contrib", "$cos", "$dmod")
    Description <- c("eigen values", "scores", "loadings", "score correlations", 
                   "distance to origin", "contribution of rows", 
                   "squared cosinus", "distance to the model")
    res1 <- cbind(Name, Description)
    cat("Principal Component Analysis by NIPALS algorithm", "\n\n")
    cat("Results available in the following objects:", "\n")
    rownames(res1) <- 1:nrow(res1)
    print(res1, print.gap=3, justify="right")
    cat("---------------------------------------------------", "\n")    
    invisible(x)
}

