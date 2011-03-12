plot.plspm.groups <-
function(x,...)
{
    test <- ifelse (x$settings[[3]]=="bootstrap", "Bootstrap test", "Permutation test")
    par(mar=c(8, 4, 4, 2))
    barplot(t(x$test[,2:3]), main=c("Path coefficients between groups",test), cex.main=1, 
            ylim=c(min(x$test[,2:3]),1.20*max(x$test[,2:3])),
            beside=T, border=NA, col=c("skyblue","blue4"), names.arg=rep("",nrow(x$test)))
    abline(h=0)
    legend("topleft", legend=colnames(x$test)[2], col="skyblue", 
           text.col="skyblue", pch=15, bty="n", pt.cex=1.5)
    legend("topright", legend=colnames(x$test)[3], col="blue4", 
           text.col="blue4", pch=15, bty="n", pt.cex=1.5)
    mtext(rownames(x$test), side=1, line=1, las=2, at=1+seq(1,3*nrow(x$test),3))
}

