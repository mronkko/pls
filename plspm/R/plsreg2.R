plsreg2 <-
function(X, Y, nc=2)
{
    X <- as.matrix(X)
    Y <- as.matrix(Y)
    if (any(is.na(X))) stop("No missing data are allowed")
    if (any(is.na(Y))) stop("No missing data are allowed")
    if (nrow(X)!=nrow(Y)) stop("different number of rows in X and Y")
    n <- nrow(X)
    p = ncol(X)
    q = ncol(Y)
    if (p<2 || q<2) stop("X and Y must have more than one column")
    if (mode(nc)!="numeric" || length(nc)!=1 || 
        nc<=1 || (nc%%1)!=0 || nc>min(n,p))
        nc <- min(n,p)   
    if (nc==n) nc<-n-1
    X.old <- scale(X)
    Y.old <- scale(Y)
    if (is.null(colnames(X)))
        colnames(X) <- paste(rep("X",p),1:p,sep="")
    if (is.null(rownames(X)))
        rownames(X) <- 1:n
    if (is.null(colnames(Y)))
        colnames(Y) <- paste(rep("Y",q),1:q,sep="")
    if (is.null(rownames(Y)))
        rownames(Y) <- 1:n
    Wh = matrix(0, p, nc)
    Uh = matrix(0, n, nc)
    Th = matrix(0, n, nc)
    Ch = matrix(0, q, nc)
    Ph = matrix(0, p, nc)
    bh <- rep(0, nc)
    RSS <- rbind(rep(n-1,q), matrix(NA, nc, q))
    PRESS <- matrix(NA, nc, q)
    Q2 <- matrix(NA, nc, q)

    for (h in 1:nc)
    {
        u.new <- Y.old[,1]# first column of Y.old
        w.old <- rep(1,p)
        iter <- 1
        repeat
        {
            w.new <- t(X.old) %*% u.new / sum(u.new^2)
            w.new <- w.new / sqrt(sum(w.new^2))# normalize w.old
            t.new <- X.old %*% w.new
            c.new <- t(Y.old) %*% t.new / sum(t.new^2)
            u.new <- Y.old %*% c.new / sum(c.new^2)
            w.dif <- w.new - w.old
            w.old <- w.new
            if (sum(w.dif^2)<1e-06 || iter==100) break
            iter <- iter + 1
        } 
        p.new <- t(X.old) %*% t.new / sum(t.new^2)

        # cross validation "leave-one-out"
        RSS[h+1,] =  colSums((Y.old - t.new%*%t(c.new))^2)
        press = matrix(0,n,q)
        for (i in 1:n)
        {
             uh.si <- Y.old[-i,1]
             wh.siold <- rep(1,p)
             itcv <- 1
             repeat
             {
                 wh.si <- t(X.old[-i,]) %*% uh.si / sum(uh.si^2)
                 wh.si <- wh.si / sqrt(sum(wh.si^2))
                 th.si <- X.old[-i,] %*% wh.si
                 ch.si <- t(Y.old[-i,]) %*% th.si / sum(th.si^2)
                 uh.si <- Y.old[-i,] %*% ch.si / sum(ch.si^2)
                 wsi.dif <- wh.si - wh.siold
                 wh.siold <- wh.si
                 if (sum(wsi.dif^2)<1e-06 || itcv==100) break
                 itcv <- itcv + 1
             }
             Yhat.si = (X.old[i,] %*% wh.si) %*% t(ch.si) 
             press[i,] = (Y.old[i,] - Yhat.si)^2
        }
        PRESS[h,] = colSums(press)
        Q2[h,] = 1 - PRESS[h,]/RSS[h,]
 
        X.old <- X.old - (t.new %*% t(p.new))
        Y.old <- Y.old - (t.new %*% t(c.new))
        Wh[,h] <- w.new
        Uh[,h] <- u.new
        Th[,h] <- t.new
        Ch[,h] <- c.new
        Ph[,h] <- p.new
        bh[h] <- t(u.new) %*% t.new 
    }
    Ws <- Wh %*% solve(t(Ph)%*%Wh)
    Bs <- Ws %*% t(Ch)# standardized regression coefficients
    Br <- diag(1/apply(X,2,sd)) %*% Bs %*% diag(apply(Y,2,sd))
    cte <- as.vector(round(apply(Y,2,mean) - apply(X,2,mean)%*%Br, 4))
    Y.hat <- X %*% Br + matrix(rep(cte,each=n),n,q)
    resids <- round(Y - Y.hat, 4)
    cor.tx <- round(cor(X,Th), 4)
    cor.ty <- round(cor(Y,Th), 4)
    R2x <- round(cor(X, Th)^2, 4)  # R2 coefficients
    R2y <- round(cor(Y, Th)^2, 4)  # R2 coefficients
    Rdx <- apply(R2x, 2, mean)
    Rdy <- apply(R2y, 2, mean)
    EV <- round(cbind(Rdx, cumsum(Rdx), Rdy, cumsum(Rdy)), 4)
    Rd.mat <- matrix(0, nc, nc)
    for (j in 1:nc)
        Rd.mat[1:j,j] <- Rdy[1:j]
    VIP <- sqrt((Wh^2) %*% Rd.mat %*% diag(p/cumsum(Rdy), nc, nc))
    Q2G <- round(1 - rowSums(PRESS)/rowSums(RSS[-nc,]), 4)
    Q2cum <- Q2
    Q2cum[1,] <- round(1-PRESS[1,]/RSS[1,], 4)
    for (i in 2:nrow(Q2))
        Q2cum[i,] <- round(1-apply(PRESS[1:i,]/RSS[1:i,], 2, prod), 4)
    Q2Gcum <- Q2G
    for (i in 1:nc)
        Q2Gcum[i] <- round(1 - prod((rowSums(PRESS)/rowSums(RSS[-nc,]))[1:i]),4)
    Q2T <- cbind(round(Q2,4),Q2G)
    Q2TC <- cbind(Q2cum,Q2Gcum)
    dimnames(Wh) <- list(colnames(X), paste(rep("w",nc),1:nc,sep=""))
    dimnames(Ws) <- list(colnames(X), paste(rep("w*",nc),1:nc,sep=""))
    dimnames(Uh) <- list(rownames(Y), paste(rep("u",nc),1:nc,sep=""))
    dimnames(Th) <- list(rownames(X), paste(rep("t",nc),1:nc,sep=""))
    dimnames(Ch) <- list(colnames(Y), paste(rep("c",nc),1:nc,sep=""))
    dimnames(Ph) <- list(colnames(X), paste(rep("p",nc),1:nc,sep=""))
    dimnames(Bs) <- list(colnames(X), colnames(Y))
    dimnames(Br) <- list(colnames(X), colnames(Y))
    dimnames(cor.tx) <- list(colnames(X), paste(rep("t",nc),1:nc,sep=""))
    dimnames(cor.ty) <- list(colnames(Y), paste(rep("t",nc),1:nc,sep=""))
    dimnames(EV) <- list(paste(rep("t",nc),1:nc,sep=""), c("R2X","R2Xcum","R2Y","R2Ycum"))
    dimnames(Y.hat) <- list(rownames(Y), colnames(Y))
    dimnames(resids) <- list(rownames(Y), colnames(Y))
    dimnames(VIP) <- list(colnames(X), paste(rep("t",nc),1:nc,sep=""))
    q2 <- c(paste(rep("Q2",q),colnames(Y),sep="."),"Q2")
    q2c <- c(paste(rep("Q2cum",q),colnames(Y),sep="."),"Q2cum")
    dimnames(Q2T) <- list(paste(rep("t",nc),1:nc,sep=""), q2)
    dimnames(Q2TC) <- list(paste(rep("t",nc),1:nc,sep=""), q2c)
    coeffs <- round(rbind(Br, INTERCEPT=cte), 4)
    res = list(x.scores=round(Th,4), x.loads=round(Ph,4), y.scores=round(Uh,4), y.loads=round(Ch,4),
              raw.wgs=round(Wh,4), weights=round(Ws,4), cor.tx=cor.tx, cor.ty=cor.ty, 
              std.coef=round(Bs,4), coeffs=coeffs, y.pred=round(Y.hat,4), resid=resids,
              expvar=EV, Q2=Q2T, Q2cum=Q2TC, VIP=round(VIP,4))
    class(res) <- "plsreg2"
    return(res)
}

