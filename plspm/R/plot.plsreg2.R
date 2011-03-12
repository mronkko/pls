plot.plsreg2 <-
function(x, ...)
{
    # Circle of correlations
    t = seq(0,2*pi,l=100)
    plot(cos(t), sin(t), type="l", main=expression(bold("Circle of Correlations on  ") * bold(list(t[1],t[2]))), 
         ylim=c(-1,1), xlab=expression("X-component  " * t[1]), ylab=expression("X-component  " * t[2]), 
         cex.main=1, xlim=c(-1,1), cex.axis=.8, col="grey")
    abline(h=seq(-1,1,.25), v=seq(-1,1,.25), col="grey", lty=3)
    abline(h=0, v=0, col="grey", lwd=2)
    points(x$cor.tx[,1], x$cor.tx[,2], pch=20, col=rep("blue",nrow(x$x.loads)))
    text(x$cor.tx[,1], x$cor.tx[,2], labels=abbreviate(rownames(x$cor.tx),7), pos=2, 
         col=rep("blue",nrow(x$x.loads)), cex=.8)
    points(x$cor.ty[,1], x$cor.ty[,2], pch=17, cex=.8, col=rep("red",nrow(x$y.loads)))
    text(x$cor.ty[,1], x$cor.ty[,2], labels=abbreviate(rownames(x$cor.ty),7), pos=2, 
         col=rep("red",nrow(x$y.loads)), cex=.8)

    # WC plot
    dev.new()
    xmi <- 1.15*min(x$weights[,1], x$y.loads[,1])
    xma <- 1.15*max(x$weights[,1], x$y.loads[,1])
    ymi <- 1.15*min(x$weights[,2], x$y.loads[,2])
    yma <- 1.15*max(x$weights[,2], x$y.loads[,2])
    plot(x$weights[,1], x$weights[,2], main=expression(bold("Loadings on weights  ") * bold(list(tilde(w),c))), 
         cex.main=1, col="blue", xlab=expression(list(tilde(w[1]),c[1])), ylab=expression(list(tilde(w[2]),c[2])), 
         pch=20, xlim=c(xmi,xma), ylim=c(ymi,yma), cex.axis=.8)
    abline(h=0, v=0, col="grey")
    text(x$weights[,1], x$weights[,2], labels=rownames(x$x.loads), cex=.7, col="blue", pos=2)
    points(x$y.loads[,1], x$y.loads[,2], col="red", pch=17, cex=.8)
    text(x$y.loads[,1], x$y.loads[,2], labels=rownames(x$y.loads), cex=.7, col="red", pos=3)

    # TT and UU plots
    dev.new()
    par(mfrow=c(1,2))
    # t1,t2 plot
    plot(x$x.scores[,1], x$x.scores[,2], main=expression(bold("Scores on components  ")* bold(list(t[1],t[2]))), 
         xlab=expression("X-component  " * t[1]), ylab=expression("X-component  " * t[2]), 
         cex.main=1, cex.axis=.8, pch=20, col="blue", cex=.8)
    abline(h=0, v=0, col="grey")
    text(x$x.scores[,1], x$x.scores[,2], labels=abbreviate(rownames(x$x.scores)), 
         cex=.7, col="blue", pos=2)
    # u1,u2 plot
    plot(x$y.scores[,1], x$y.scores[,2], main=expression(bold("Scores on components  ")* bold(list(u[1],u[2]))), 
         xlab=expression("Y-component  " * u[1]), ylab=expression("Y-component  " * u[2]), 
         cex.main=1, cex.axis=.8, pch=17, cex=.7, col="red")
    abline(h=0, v=0, col="grey")
    text(x$y.scores[,1], x$y.scores[,2], labels=abbreviate(rownames(x$y.scores)), 
         cex=.7, col="red", pos=2)

    # TU plots
    dev.new()
    par(mfrow=c(1,2))
    # t1,u1 plot
    plot(x$x.scores[,1], x$y.scores[,1], main=expression(bold("Scores on components  ")* bold(list(t[1],u[1]))), 
         xlab=expression("X-component  " * t[1]), ylab=expression("Y-component  " * u[1]), 
         cex.main=1, cex.axis=.8, type="n")
    abline(h=0, v=0, col="grey")
    points(x$x.scores[,1], x$y.scores[,1], col="green3", pch=20, cex=.8)
    text(x$x.scores[,1], x$y.scores[,1], labels=abbreviate(rownames(x$y.scores)), 
         cex=.7, col="green3", pos=2)
    # t2,u2 plot
    plot(x$x.scores[,2], x$y.scores[,2], main=expression(bold("Scores on components  ")* bold(list(t[2],u[2]))), 
         xlab=expression("X-component  " * t[2]), ylab=expression("Y-component  " * u[2]), 
         cex.main=1, cex.axis=.8, type="n")
    abline(h=0, v=0, col="grey")
    points(x$x.scores[,2], x$y.scores[,2], col="green3", pch=20, cex=.8)
    text(x$x.scores[,2], x$y.scores[,2], labels=abbreviate(rownames(x$y.scores)),
         cex=.7, col="green3", pos=2)
}

