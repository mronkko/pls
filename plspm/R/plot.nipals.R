plot.nipals <-
function(x, ...)
{
    # screeplot
    par(mfrow=c(1,2))
    barplot(x$values[,1], main="Screeplot of eigenvalues", cex.main=1, 
            ylab="values", border=NA, names.arg=rownames(x$values))
    barplot(x$values[,2], main="Screeplot of eigenvalues", cex.main=1, 
            border=NA, names.arg=rownames(x$values), ylim=c(0,100),
            ylab="Percentage of explained variance")   
    points(x$values[,3], pch=20, col="blue")
    lines(x$values[,3], col="blue")
    # Circle of correlations
    dev.new()
    t = seq(0,2*pi,l=100)
    plot(cos(t), sin(t), type="l", main="Circle of Correlations", 
         xlab=expression("Component  " * t[1]), ylab=expression("Component  " * t[2]), 
         cex.main=1, xlim=c(-1,1), ylim=c(-1,1), col="grey", cex.axis=.8)
    abline(h=seq(-1,1,.25), v=seq(-1,1,.25), col="grey", lty=3)
    abline(h=0, v=0, col="grey")
    points(x$cor.sco[,1], x$cor.sco[,2], pch=17, col="blue", cex=.8)
    text(x$cor.sco[,1], x$cor.sco[,2], labels=abbreviate(rownames(x$loadings)),  
         pos=2, col="blue", cex=.7)
    # Graphic of components TT
    dev.new()
    par(mfrow=c(1,2))
    plot(x$scores[,1], x$scores[,2], main=expression(bold("Graphic of components  " * t[1] *","* t[2])), 
         xlab=expression("Component  " * t[1]), ylab=expression("Component  " * t[2]),
         cex.main=1, cex.axis=.8, pch=20, col="blue")
    abline(h=0, v=0, col="grey")
    text(x$scores[,1], x$scores[,2], labels=abbreviate(rownames(x$scores)), col="blue", pos=2, 
         offset=.2, cex=.7)
    # Graphic of loadings PP
    plot(x$loadings[,1], x$loadings[,2], main=expression(bold("Graphic of loadings  " * p[1] *","* p[2])), 
         xlab=expression("Loading  " * p[1]), ylab=expression("Loading  " * p[2]),
         cex.main=1, cex.axis=.8, pch=20, col="red")
    abline(h=0, v=0, col="grey")
    text(x$loadings[,1], x$loadings[,2], labels=abbreviate(rownames(x$loadings)), col="red", pos=2, 
         offset=.2, cex=.7)
}

