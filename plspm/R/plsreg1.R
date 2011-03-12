`plsreg1` <-
function(x, y, nc=2, cv=FALSE)
{
    # ============ checking arguments ============
    X <- as.matrix(x)
    n <- nrow(X)
    p <- ncol(X)
    if (p < 2) stop("'x' must contain more than one column")
    if (is.null(colnames(X)))
        colnames(X) <- paste(rep("X",p),1:p,sep="")
    if (is.null(rownames(X)))
        rownames(X) <- 1:n
    if (!is.vector(y))
        stop("Invalid argument 'y'. Must be a numeric vector")
    if (any(is.na(y))) 
        stop("'y' must not contain missing values")    
    Y <- as.matrix(y)
    if (nrow(X)!=nrow(Y)) 
        stop("x and y have different number of rows")
    if (is.null(colnames(Y)))
        colnames(Y) <- "Y"
    if (is.null(rownames(Y)))
        rownames(Y) <- 1:n
    if (mode(nc)!="numeric" || length(nc)!=1 || 
        nc<=1 || (nc%%1)!=0 || nc>min(n,p))
        nc <- min(n,p)   
    if (nc==n) nc<-n-1
    if (!is.logical(cv)) cv<-FALSE
    if (any(is.na(X))) na.miss<-TRUE else na.miss<-FALSE       
    if (na.miss) cv <- FALSE
    if (cv) nc <- min(n,p) 
    # ============ setting inputs ==============
    Xx <- scale(X) 
    Yy <- scale(Y)
    X.old <- Xx
    Y.old <- Yy
    Th <- matrix(NA, n, nc)# matrix of X-scores
    Ph <- matrix(NA, p, nc)# matrix of X-loadings
    Wh <- matrix(NA, p, nc)# matrix of raw-weights
    Uh <- matrix(NA, n, nc)# matrix of Y-scores
    ch <- rep(NA, nc)# vector of y-loadings
    RSS <- c(n-1, rep(NA, nc))
    PRESS <- rep(NA, nc)
    Q2 <- rep(NA, nc)
    Hot <- matrix(NA, n, nc)# matrix of T2 Hotelling  
    hlim <- rep(NA, nc)# matrix of T2 thresholds
    # ============ pls regression algorithm ==============
    w.old <- rep(1,p)
    t.new <- rep(1,n)
    p.new <- rep(NA, p)
    h <- 1
    repeat
    {
        if (na.miss)# missing data
        {
            for (j in 1:p) {
                i.exist <- which(complete.cases(X[,j]))
                w.old[j] <- sum(X.old[i.exist,j] * Y.old[i.exist])
            }
            w.new <- w.old / sqrt(sum(w.old^2))       
            for (i in 1:n) {
                j.exist <- which(complete.cases(X[i,]))
                t.new[i] <- sum(X.old[i,j.exist] * w.new[j.exist])
            }
            for (j in 1:p) {
                i.exist <- intersect(which(complete.cases(X[,j])), which(complete.cases(t.new)))
                p.new[j] <- sum(X.old[i.exist,j] * t.new[i.exist]) / sum(t.new[i.exist]^2)
            }
        }
        if (!na.miss) # no missing data
        {
            w.old <- t(X.old) %*% Y.old / sum(Y.old^2)
            w.new <- w.old / sqrt(sum(w.old^2)) # normalization
            t.new <- X.old %*% w.new
            p.new <- t(X.old) %*% t.new / sum(t.new^2)
        }
        c.new <- t(Y.old) %*% t.new / sum(t.new^2)
        u.new <- Y.old / as.vector(c.new)
        if (!na.miss) 
        {
            # cross validation "leave-one-out"
            RSS[h+1] <-  sum((Y.old - t.new%*%c.new)^2)
            press <- rep(0,n)
            for (i in 1:n)
            {
                 Xy.aux <- t(X.old[-i,]) %*% Y.old[-i]
                 wh.si <- Xy.aux %*% sqrt(solve(t(Xy.aux)%*%Xy.aux))
                 th.si <- X.old[-i,] %*% wh.si
                 ch.si <- t(Y.old[-i]) %*% th.si %*% solve(t(th.si)%*%th.si)
                 ch.si <- as.vector(ch.si)
                 Yhat.si <- ch.si * X.old[i,] %*% wh.si
                 press[i] <- (Y.old[i] - Yhat.si)^2
            }
            PRESS[h] = sum(press)
            Q2[h] = 1 - PRESS[h]/RSS[h]
        }
        Y.old <- Y.old - t.new%*%c.new# deflate y.old
        X.old <- X.old - (t.new %*% t(p.new))# deflate X.old
        Th[,h] <- t.new
        Ph[,h] <- p.new
        Wh[,h] <- w.new
        Uh[,h] <- u.new
        ch[h] <- c.new   
        Hot[,h] <- (n/(n-1)) * t.new^2 / (sum(t.new^2)/(n-1))
        hlim[h] <- qf(.95,h,n-h) * (h*(n^2-1))/(n*(n-h)) 
        if (!na.miss)
           if (cv) 
              if (Q2[h]<0.0975 || h==nc) break
        if (!cv)
           if (h==nc) break
        h <- h + 1
    }
    Th <- round(Th[,1:h], 4)
    Ph <- round(Ph[,1:h], 4)
    Wh <- round(Wh[,1:h], 4)
    Uh <- round(Uh[,1:h], 4)
    ch <- round(ch[1:h], 4)
    Ws <- round(Wh[,1:h] %*% solve(t(Ph[,1:h])%*%Wh[,1:h]), 4)# modified weights
    Bs <- round(as.vector(Ws %*% ch[1:h]), 4) # std beta coeffs    
    if (!na.miss) 
    {   
        Br <- round(Bs * (rep(sd(Y),p)/apply(X,2,sd)), 4)   # beta coeffs
        cte <- as.vector(round(mean(y) - Br%*%apply(X,2,mean), 4))# intercept
        y.hat <- round(as.vector(X%*%Br + cte), 4)# y predicted
        q2cum <- rep(NA,h)
        for (k in 1:h)
            q2cum[k] <- prod(PRESS[1:k]) / prod(RSS[1:k])
        Q2cum <- round((1 - q2cum), 4)
        Q2cv <- round(cbind(PRESS[1:h], RSS[1:h], Q2[1:h], rep(0.0975,h), Q2cum), 4)
        dimnames(Q2cv) <- list(1:h, c("PRESS","RSS","Q2","LimQ2","Q2cum"))
        cor.sco<-round(cor(cbind(Xx, y=Yy), Th), 4)
    } else
    {
        mu.x <- attributes(Xx)$'scaled:center'
        sd.x <- attributes(Xx)$'scaled:scale'
        X.hat <- Th%*%t(Ph) %*% diag(sd.x,p,p) + matrix(rep(mu.x,each=n),n,p)
        Br <- round(Bs * (rep(sd(Y),p)/sd.x), 4)   # beta coeffs
        cte <- as.vector(round(mean(y) - Br%*%mu.x, 4))# intercept
        y.hat <- round(as.vector(X.hat%*%Br + cte), 4)
        cor.sco <- matrix(NA, p+1, nc)
        for (j in 1:p) {
            i.exist <- which(complete.cases(X[,j]))
            cor.sco[j,] <- round(cor(Xx[i.exist,j], Th[i.exist,]), 4)
        }
        cor.sco[p+1,] <- round(cor(Yy, Th), 4)
    }
    resid <- round(as.vector(Y - y.hat), 4)# residuals
    R2 <- round(as.vector(cor(Th, Yy))^2, 4)  # R2 coefficients
    T2hot <- round(rbind(hlim[1:h], t(apply(Hot[,1:h],1,cumsum))), 4)# Hotelling T2 

    dimnames(Wh) <- list(colnames(X), paste(rep("w",h),1:h,sep=""))
    dimnames(Ws) <- list(colnames(X), paste(rep("w*",h),1:h,sep=""))    
    dimnames(Th) <- list(rownames(X), paste(rep("t",h),1:h,sep=""))
    dimnames(Ph) <- list(colnames(X), paste(rep("p",h),1:h,sep=""))
    dimnames(Uh) <- list(rownames(Y), paste(rep("u",h),1:h,sep=""))
    names(ch) <- paste(rep("c",h),1:h,sep="")
    dimnames(T2hot) <- list(c("T2",rownames(X)), paste(rep("H",h),1:h,sep=""))
    names(Bs) <- colnames(X)
    names(Br) <- colnames(X)
    names(resid) <- rownames(Y)
    names(y.hat) <- rownames(Y)
    names(R2) <- paste(rep("t",h),1:h,sep="")
    dimnames(cor.sco) <- list(c(colnames(X), colnames(Y)), colnames(Th))
    if (!na.miss) 
        res <- list(scores=Th, x.loads=Ph, y.loads=ch, u.scores=Uh, raw.wgs=Wh,
                    mod.wgs=Ws, std.coef=Bs, coeffs=c(Intercept=cte,Br), R2=R2, 
                    y.pred=y.hat, resid=resid, cor.sco=cor.sco, T2=T2hot, Q2=Q2cv, y=Y)    
    if (na.miss)
        res <- list(scores=Th, x.loads=Ph, y.loads=ch, u.scores=Uh, raw.wgs=Wh,
                    mod.wgs=Ws, std.coef=Bs, coeffs=c(Intercept=cte,Br), R2=R2, 
                    y.pred=y.hat, resid=resid, cor.sco=cor.sco, T2=T2hot, y=Y)         
    class(res) <- "plsreg1"
    return(res)
}

