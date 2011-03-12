`nipals` <-
function(x, nc=2, scaled=TRUE)
{
    # checking arguments
    X <- as.matrix(x)
    n <- nrow(X)
    p <- ncol(X)
    if (!n || !p) stop("dimension 0 in 'x'")
    if (p==1) stop("x must be a numeric matrix or data frame")
    if (is.null(colnames(X)))
        colnames(X) <- paste(rep("X",p),1:p,sep="")
    if (is.null(rownames(X)))
        rownames(X) <- rep(1:n)
    if (!is.logical(scaled)) scaled <- TRUE
    if (scaled) X<-scale(X) else X<-scale(X,scale=FALSE)    
    if (mode(nc)!="numeric" || length(nc)!=1 || 
        nc<=1 || (nc%%1)!=0 || nc>min(n,p))
        nc <- 2
    if (nc==n) nc <- n-1
    if (any(is.na(X))) na.miss<-TRUE else na.miss<-FALSE       
    # setting inputs
    X.old = X    
    Th = matrix(NA, n, nc)# scores
    Ph = matrix(NA, p, nc)# loadings
    eig.vals = rep(NA,nc)# eigenvalues
    for (h in 1:nc)
    {
        th.new <- X.old[,1]
        ph.old <- rep(1,p)
        ph.new <- rep(1,p)
        iter <- 1
        repeat
        {
            if (na.miss)# missing data
            {
                for (j in 1:p)
                {
                    i.exist <- intersect(which(complete.cases(X[,j])), which(complete.cases(th.new)))
                    ph.new[j] <- sum(X.old[i.exist,j] * th.new[i.exist]) / sum(th.new[i.exist]^2)
                }
                ph.new <- ph.new / sqrt(sum(ph.new^2))       
                for (i in 1:n)
                {
                    j.exist <- which(complete.cases(X[i,]))
                    th.new[i] <- sum(X.old[i,j.exist] * ph.new[j.exist]) / sum(ph.new[j.exist]^2)
                }
            } else
            {
                ph.new <- t(X.old) %*% th.new / sum(th.new^2)
                ph.new <- ph.new / sqrt(sum(ph.new^2))
                th.new <- X.old %*% ph.new
            }
            ph.aux <- sum((ph.new - ph.old)^2)
            ph.old <- ph.new
            if (ph.aux < 1e-06 || iter==100) break
               iter <- iter + 1
        }
        Th[,h] <- round(th.new, 4)
        Ph[,h] <- round(ph.new, 4)
        X.new <- X.old - th.new %*% t(ph.new)
        X.old <- X.new
        eig.vals[h] <- sum(th.new^2)/(n-1)
    }
    # eigenvalues
    eig.perc <- round(100*eig.vals/p, 2)
    eigs <- data.frame(values=round(eig.vals,4), val.perc=eig.perc, val.cum=cumsum(eig.perc))
    rownames(eigs) <- paste(rep("v",nc),1:nc,sep="")
    dimnames(Th) <- list(rownames(X), paste(rep("t",nc),1:nc,sep=""))
    dimnames(Ph) <- list(colnames(X), paste(rep("p",nc),1:nc,sep=""))
    # interpretation tools
    if (na.miss) {
        cor.sco <- matrix(NA, p, nc)
        for (j in 1:p) {
            i.exist <- which(complete.cases(X[,j]))
            cor.sco[j,] <- round(cor(X[i.exist,j], Th[i.exist,]), 4)
        }
        dimnames(cor.sco) <- list(colnames(X), colnames(Th))
    } else
        cor.sco <- round(cor(X, Th), 4)
    ConInd <- round((100/n)*Th^2 %*% diag(1/eig.vals),4)# individuals contribution
    dimnames(ConInd) <- list(rownames(X), paste(rep("ctr",nc),1:nc,sep=".") )
    D.proj <- Th^2# square distance of projections
    D.orig <- rowSums(X^2, na.rm=T)# square distance of centered data
    Cos.inds <- matrix(0,n,nc)# square cosinus
    for (i in 1:n)
        Cos.inds[i,] <- round(D.proj[i,]/D.orig[i], 4)
    dimnames(Cos.inds) <- list(rownames(X), paste(rep("cos2",nc),1:nc,sep=".") )
    # Distance to the model: DModX 
    DMX2 <- matrix(NA, n, nc)
    DMX2[,1] <- (D.orig * (1-Cos.inds[,1])) / (p-eig.vals[1])
    for (j in 2:nc)
         DMX2[,j] <- (D.orig * (1- rowSums(Cos.inds[,1:j]))) / (p-sum(eig.vals[1:j]))
    DMX = round(DMX2, 4)# Distance to the Model X (normalized)
    dimnames(DMX) <- list(rownames(X), colnames(Th))
    # results 
    res <- list(values=eigs, scores=Th, loadings=Ph, cor.sco=cor.sco, 
           disto=round(D.orig,4), contrib=ConInd, cos=Cos.inds, dmod=DMX)
    class(res) <- "nipals"
    return(res)
}

