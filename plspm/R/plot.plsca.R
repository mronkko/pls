plot.plsca <-
function(x, ...)
{
    # Circle of correlations
    par (mfrow=c(1,2))
    t = seq(0,2*pi,l=100)
    plot(cos(t), sin(t), type="l", main="Circle of Correlations (X-components)", 
          xlab=expression("X-component  " * t[1]), ylab=expression("X-component  " * t[2]), 
         cex.main=1, xlim=c(-1,1), ylim=c(-1,1), col="grey")
    abline(h=seq(-1,1,.25), v=seq(-1,1,.25), col="grey", lty=3)
    abline(h=0, v=0, col="grey")
    points(x$cor.xt[,1], x$cor.xt[,2], pch=20, col="blue")
    points(x$cor.yt[,1], x$cor.yt[,2], pch=17, col="red", cex=.8)
    text(x$cor.xt[,1], x$cor.xt[,2], labels=rownames(x$cor.xt), pos=3, col="blue", cex=.7)
    text(x$cor.yt[,1], x$cor.yt[,2], labels=rownames(x$cor.yt), pos=3, col="red", cex=.7)
    plot(cos(t), sin(t), type="l", main="Circle of Correlations (Y-components)",
          xlab=expression("Y-component  " * u[1]), ylab=expression("Y-component  " * u[2]), 
         cex.main=1, xlim=c(-1,1), ylim=c(-1,1), col="grey")
    abline(h=seq(-1,1,.25), v=seq(-1,1,.25), col="grey", lty=3)
    abline(h=0, v=0, col="grey")
    points(x$cor.yu[,1], x$cor.yu[,2], pch=17, col="red", cex=.8)
    points(x$cor.xu[,1], x$cor.xu[,2], pch=20, col="blue")
    text(x$cor.yu[,1], x$cor.yu[,2], labels=rownames(x$cor.yu), pos=3, col="red", cex=.7)
    text(x$cor.xu[,1], x$cor.xu[,2], labels=rownames(x$cor.xu), pos=3, col="blue", cex=.7)

    # plot of scores
    dev.new()
    par(mfrow=c(1,2))
    plot(x$x.scores[,1], x$x.scores[,2], main=expression(bold("Graphic of PLS components  " * t[1] *","* t[2])), 
         xlab=expression("X-component  " * t[1]), ylab=expression("X-component  " * t[2]),
         cex.main=1, cex.axis=.8, pch=20, col="blue")
    abline(h=0, v=0, col="grey")
    text(x$x.scores[,1], x$x.scores[,2], labels=abbreviate(rownames(x$x.scores)), col="blue", pos=2, 
         offset=.2, cex=.7)
    plot(x$y.scores[,1], x$y.scores[,2], main=expression(bold("Graphic of PLS components  " * u[1] *","* u[2])), 
         xlab=expression("Y-component  " * u[1]), ylab=expression("Y-component  " * u[2]),
         cex.main=1, cex.axis=.8, pch=17, col="red", cex=.8)
    abline(h=0, v=0, col="grey")
    text(x$y.scores[,1], x$y.scores[,2], labels=abbreviate(rownames(x$y.scores)), col="red", pos=2, 
         offset=.2, cex=.7)

    # plot of R2 and communalities
    dev.new()
    par(mfrow=c(2,2), mar=c(2.5, 3, 3, 1.5))
    barplot(x$R2X, beside=T, border=NA, main="Explained variance of X-scores", cex.main=1)
    barplot(x$R2Y, beside=T, border=NA, main="Explained variance of Y-scores", cex.main=1)
    barplot(x$com.xu, beside=T, border=NA, main="Communality of X with U", cex.main=1)
    barplot(x$com.yt, beside=T, border=NA, main="Communality of Y with T", cex.main=1)

    # plot T-vs-U
    dev.new()
    (z1 <- line(cbind(x$x.scores[,1], x$y.scores[,1])))
    par(mfrow=c(1,2))
    plot(x$x.scores[,1], x$y.scores[,1], main=expression(bold("Graphic of PLS components  " * t[1] *","* u[1])), 
         xlab=expression("X-component  " * t[1]), ylab=expression("Y-component  " * u[1]),
         cex.main=1, cex.axis=.8, pch=20, col="green4")
    abline(h=0, v=0, col="grey")
    abline(coef(z1),col="green")
    text(x$x.scores[,1], x$y.scores[,1], labels=abbreviate(rownames(x$x.scores)), col="green4", pos=2, 
         offset=.2, cex=.7) 
    (z2 <- line(cbind(x$x.scores[,2], x$y.scores[,2])))
    plot(x$x.scores[,2], x$y.scores[,2], main=expression(bold("Graphic of PLS components  " * t[2] *","* u[2])), 
         xlab=expression("X-component  " * t[2]), ylab=expression("Y-component  " * u[2]),
         cex.main=1, cex.axis=.8, pch=20, col="green4")
    abline(h=0, v=0, col="grey")
    abline(coef(z2),col="green")
    text(x$x.scores[,2], x$y.scores[,2], labels=abbreviate(rownames(x$x.scores)), col="green4", pos=2, 
         offset=.2, cex=.7)   
}

