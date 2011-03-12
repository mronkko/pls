plot.plsreg1 <-
function(x, ...)
{
    # Circle of correlations
    t = seq(0,2*pi,l=100)
    plot(cos(t), sin(t), type="l", main="Circle of Correlations", 
          xlab=expression("X-component  " * t[1]), ylab=expression("X-component  " * t[2]), 
         cex.main=1, xlim=c(-1,1), ylim=c(-1,1), col="grey")
    abline(h=seq(-1,1,.25), v=seq(-1,1,.25), col="grey", lty=3)
    abline(h=0, v=0, col="grey")
    points(x$cor.sco[,1], x$cor.sco[,2], pch=20, col=c(rep("blue",nrow(x$x.loads)),"red"))
    text(x$cor.sco[,1], x$cor.sco[,2], labels=abbreviate(rownames(x$cor.sco)), pos=2, 
         col=c(rep("blue",nrow(x$x.loads)),"red"), cex=.7)
    # Graphic of components TT
    dev.new()
    par(mfrow=c(1,2))
    plot(x$scores[,1], x$scores[,2], main=expression(bold("Graphic of PLS components  " * t[1] *","* t[2])), 
         xlab=expression("X-component  " * t[1]), ylab=expression("X-component  " * t[2]),
         cex.main=1, cex.axis=.8, pch=20, col="blue")
    abline(h=0, v=0, col="grey")
    text(x$scores[,1], x$scores[,2], labels=abbreviate(rownames(x$scores)), col="blue", pos=2, 
         offset=.2, cex=.7)
    # Graphic of components TU
    plot(x$scores[,1], x$u.scores[,1], main=expression(bold("Graphic of PLS components  " * t[1] *","* u[1])), 
          xlab=expression("X-component  " * t[1]), ylab=expression("Y-component  " * u[1]),
         cex.main=1, cex.axis=.8, pch=20, col="blue",   )
    abline(h=0, v=0, col="grey")
    text(x$scores[,1], x$u.scores[,1], labels=abbreviate(rownames(x$scores)), col="blue", pos=2, 
         offset=.2, cex=.7)
    # T2 Hotelling Confidence ellipse
    n <- nrow(x$scores)
    hot.lim <- ((n-1)/n)* (2*(n^2-1))/(n*(n-2)) * qf(.95,2,n-2)
    a <- sqrt(hot.lim * sum(x$scores[,1]^2) / (n-1))
    b <- sqrt(hot.lim * sum(x$scores[,2]^2) / (n-1))
    ellipsePoints <- function(a, b, n=201)
    {
        ## a, b : length of half axes in (x,y) direction
        ## n    : number of points
        B <- min(a,b)
        A <- max(a,b)
        d2 <- (A-B)*(A+B)                   #= A^2 - B^2
        phi <- 2*pi*seq(0,1, len = n)
        r <- a*b / sqrt(B^2 + d2 * sin(phi)^2)
        xy <- r * cbind(cos(phi), sin(phi))
        xy %*% rbind(c(cos(0),sin(0)), c(-sin(0),cos(0))) + cbind(rep(0,n),rep(0,n))
    }
    ep1 = ellipsePoints(a, b)
    dev.new()
    plot(x$scores[,1], x$scores[,2], main=expression(bold(" Ellipse Hotelling  " *T^2 * " (0.05)")), 
         type="n", cex.main=1, cex.lab=.8, xlab=expression(t[1]), ylab=expression(t[2]),
          cex.axis=.8, xlim=1.10*c(-a,a), ylim=1.10*c(-b,b) ) 
    abline(h=0, v=0, col="grey")
    points(x$scores[,1], x$scores[,2], pch=20, col="blue")
    points(ep1, type="l", col="purple", lty=1)
    text(x$scores[,1], x$scores[,2], labels=abbreviate(rownames(x$scores)), col="blue", pos=2, cex=.7)
    # Comparison Y -vs- Y.hat
    dev.new()
    (z <- line(cbind(x$y, x$y.pred)))
    plot(x$y, x$y.pred, main="Y-observed -vs- Y-predicted", type="n",
         xlab="Y", ylab=expression(hat(Y)), cex.main=1, cex.axis=.8)
    abline(coef(z),col="green")
    points(x$y, x$y.pred, pch=20, col="green4")
    text(x$y, x$y.pred, labels=abbreviate(rownames(x$y)), col="green4", pos=2, cex=.7)
}

