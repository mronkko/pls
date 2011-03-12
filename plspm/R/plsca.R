`plsca` <-
function(X, Y, nc=NULL, scaled=TRUE)
{
    # ============ checking arguments ============
    X <- as.matrix(X)
    Y <- as.matrix(Y)
    if (nrow(X)!=nrow(Y))
        stop("Different number of rows in X and Y")
    n <- nrow(X)
    p <- ncol(X)
    q <- ncol(Y) 
    if (p==1) stop("X must contain more than one variable")
    if (q==1) stop("Y must contain more than one variable")
    if (any(is.na(X)) | any(is.na(Y))) stop("No missing data are allowed")
    if (is.null(colnames(X)))
        colnames(X) <- paste(rep("X",p),1:p,sep="")
    if (is.null(rownames(X)))
        rownames(X) <- rep(1:n)
    if (is.null(colnames(Y)))
        colnames(Y) <- paste(rep("Y",q),1:q,sep="")
    if (is.null(rownames(Y)))
        rownames(Y) <- rep(1:n)        
    if (!is.logical(scaled)) scaled<-TRUE
    if (!is.null(nc)) 
    {
        if (mode(nc)!="numeric" || length(nc)!=1 || 
            nc<=1 || (nc%%1)!=0 || nc>min(n,p,q))
            nc <- min(n,p,q)   
    } else
        nc <- min(n,p,q)
    # ============ setting inputs ==============    
    if (scaled) X<-scale(X) else X<-scale(X,scale=F)
    if (scaled) Y<-scale(Y) else Y<-scale(Y,scale=F)
    X.old <- scale(X)
    Y.old <- scale(Y)
    Xt <- matrix(NA,n,nc)
    Yu <- matrix(NA,n,nc)
    Ah <- matrix(NA,p,nc)
    Bh <- matrix(NA,q,nc)
    Ch <- matrix(NA,p,nc)
    Eh <- matrix(NA,q,nc)
    # ============ pls canonical analysis ==============
    for (h in 1:nc)
    {
        th <- X.old[,1]
        uh <- Y.old[,1]
        ah.old <- rep(1,p)
        iter <- 1
        repeat
        {
            ah <- t(X.old)%*%uh / sum(uh^2)
            ah <- ah / sqrt(sum(ah^2))
            th <- X.old %*% ah 
            bh <- t(Y.old)%*%th / sum(th^2)
            bh <- bh / sqrt(sum(bh^2))
            uh <- Y.old %*% bh 
            ah.dif <- sum((ah - ah.old)^2)
            if (ah.dif<1e-06 || iter==100) break
            iter <- iter + 1
            ah.old <- ah
        }        
        ch <- t(X.old) %*% th / sum(th^2)
        eh <- t(Y.old) %*% uh / sum(uh^2)
        X.old <- X.old - th%*%t(ch)
        Y.old <- Y.old - uh%*%t(eh)
        Xt[,h] <- round(th, 4)
        Yu[,h] <- round(uh, 4)
        Ah[,h] <- round(ah, 4)
        Bh[,h] <- round(bh, 4)
        Ch[,h] <- round(ch, 4)
        Eh[,h] <- round(eh, 4)
    }
    dimnames(Ch) <- list(colnames(X), paste(rep("c",nc),1:nc,sep=""))
    dimnames(Eh) <- list(colnames(Y), paste(rep("e",nc),1:nc,sep=""))
    dimnames(Xt) <- list(rownames(X), paste(rep("t",nc),1:nc,sep=""))
    dimnames(Yu) <- list(rownames(Y), paste(rep("u",nc),1:nc,sep=""))
    As <- round(Ah %*% solve(t(Ch)%*%Ah), 4)
    Bs <- round(Bh %*% solve(t(Eh)%*%Bh), 4)
    dimnames(As) <- list(colnames(X), paste(rep("t",nc),1:nc,sep=""))
    dimnames(Bs) <- list(colnames(Y), paste(rep("u",nc),1:nc,sep="")) 
    XT <- round(scale(X) %*% As, 4)
    YU <- round(scale(Y) %*% Bs, 4)
    dimnames(XT) <- list(rownames(X), paste(rep("t",nc),1:nc,sep=""))
    dimnames(YU) <- list(rownames(Y), paste(rep("u",nc),1:nc,sep=""))
    cor.xt <- round(cor(X,Xt), 3)
    cor.yu <- round(cor(Y,Yu), 3)
    cor.tu <- round(cor(Xt,Yu), 3)
    cor.xu <- round(cor(X,Yu), 3)
    cor.yt <- round(cor(Y,Xt), 3)
    R2X <- round(t(apply(cor.xt^2, 1, cumsum)), 3)
    R2Y <- round(t(apply(cor.yu^2, 1, cumsum)), 3)
    Com.XU <- round(t(apply(cor.xu^2, 1, cumsum)), 3)
    Com.YT <- round(t(apply(cor.yt^2, 1, cumsum)), 3)
    res <- list(x.scores=XT, x.wgs=As, x.loads=Ch, y.scores=YU, y.wgs=Bs, y.loads=Eh,
                cor.xt=cor.xt, cor.yu=cor.yu, cor.tu=cor.tu, cor.xu=cor.xu, cor.yt=cor.yt,
                R2X=R2X, R2Y=R2Y, com.xu=Com.XU, com.yt=Com.YT)
    class(res) <- "plsca"
    return(res)
}

