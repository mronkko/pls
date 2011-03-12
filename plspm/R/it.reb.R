it.reb <-
function(pls, hclus.res, nk, Y=NULL, stop.crit=0.005, iter.max=100)
{
    # ========================== it.reb function ==========================
    # Performs iterative steps of Response-Based Unit Segmentation  
    # in Partial Least Squares Path Modeling (PLS-PM)        
    # =========================== ARGUMENTS ==============================
    # pls: object of class "plspm"
    # hclus.res: object of class "hclust" obtained from function "res.clus"
    # nk: an integer larger than 1 indicating the number of classes 
    # Y: optional data matrix used when pls$data is null
    # stop.crit: stop criterion number (must be between 0 and 1)
    # iter.max: maximum number of iterations (must be an integer)

    # ==================== Checking function arguments ===================
    if (class(pls)!="plspm") 
        stop("argument 'pls' must be an object of class 'plspm'")
    if (any(pls$model$modes!="A"))# checking reflective modes
        stop("REBUS only works for reflective modes")
    if (!pls$model$scaled)# checking scaled data
        stop("REBUS only works with scaled='TRUE'")
    if (missing(hclus.res))
        stop("argument 'hclus.res' is missing")
    if (class(hclus.res)!="hclust")
        stop("argument 'hclus.res' must be an object of class 'hclust'")
    if (missing(nk))
        stop("argument 'nk' (number of classes) is missing")
    if (mode(nk)!="numeric" || length(nk)!=1 || 
        nk<=1 || (nk%%1)!=0)
        stop("Invalid number of classes 'nk'. Must be an integer larger than 1")
    if (!is.null(Y)) # if Y available
    {
        if (is.null(pls$data))
        {
            if (!is.matrix(Y) && !is.data.frame(Y))
                stop("Invalid object 'Y'. Must be a numeric matrix or data frame.")
            if (nrow(Y)!=nrow(pls$latents))
                stop("Argument 'pls' and 'Y' are incompatible. Different number of rows.")
        }
    } else { # if no Y
        if (is.null(pls$data)) 
            stop("Argument 'Y' is missing. No dataset available.")
    }
    if (mode(stop.crit)!="numeric" || length(stop.crit)!=1 || 
        stop.crit<0 || stop.crit>=1)
    {
        warning("Invalid stop criterion 'stop.crit'. Deafult value 0.005 is used")
        stop.crit <- 0.005
    }
    if (mode(iter.max)!="numeric" || length(iter.max)!=1 || 
        iter.max<=1 || (iter.max%%1)!=0)
    {
        warning("Invalid number of maximum iterations 'iter.max'. Deafult value 100 is used")
        iter.max <- 100
    }

    # ========================== INPUTS SETTING ==========================
    IDM <- pls$model$IDM# Inner Design Matrix
    blocks <- pls$model$blocks# cardinality of blocks
    scheme <- pls$model$scheme# inner weighting scheme
    modes <- pls$model$modes# measurement modes
    scaled <- pls$model$scaled# type of scaling
    plsr <- pls$model$plsr# pls-regression
    tol <- pls$model$tol# tolerance criterion
    iter <- pls$model$iter# max num iterations
    outer <- pls$model$outer
    if (plsr) 
        warning("path coefficients will be calculated with OLS regression")
    plsr <- FALSE
    blocklist <- outer
    for (k in 1:length(blocks))
         blocklist[[k]] <- rep(k,blocks[k])
    blocklist <- unlist(blocklist)
    # data matrix DM
    if (!is.null(pls$data)) {
        DM <- pls$data
    } else {         
        # building data matrix 'DM'
        DM <- matrix(NA, nrow(pls$latents), sum(blocks))
        for (k in 1:nrow(IDM))
            DM[,which(blocklist==k)] <- as.matrix(Y[,outer[[k]]])
        dimnames(DM) <- list(rownames(pls$latents), names(pls$out.weights))
    }
    lvs <- nrow(IDM)# number of LVs
    lvs.names <- rownames(IDM)# names of LVs
    mvs <- sum(blocks)# number of MVs
    mvs.names <- colnames(DM)
    N <- nrow(DM)# number of observations
    endo <- rowSums(IDM)
    endo[endo!=0] <- 1  # indicator of endog LVs
    n.end <- sum(endo)# number of enfod LVs
    # data scaling
    sd.X <- sqrt((N-1)/N) * apply(DM, 2, sd)
    X <- scale(DM, scale=sd.X) 
    # initial partition
    ini.part <- cutree(hclus.res, nk)# cutting dendrogram in 'nk' clusters
    nclus <- nlevels(factor(ini.part))  # number of clusters

    # =============== creating objects to store results ======================
    w.locals <- as.list(1:nclus)# outer.weights
    Y.locals <- as.list(1:nclus)# std latent variables 
    loads.locals <- as.list(1:nclus)# loadings
    path.locals <- as.list(1:nclus)# path coefficients
    R2.locals <- as.list(1:nclus)# R2
    comu.locals <- as.list(1:nclus)# mvs communalities
    outres.locals <- as.list(1:nclus)# communality residuals
    innres.locals <- as.list(1:nclus)# structural residuals
    CM.locals <- matrix(NA,N,nclus)# CM distance of each unit from each model

    # ====================== iterative process =====================
    old.clas <- ini.part
    iter.ch <- 0
    repeat 
    {
        # define MV matrix for each initial class
        split.DM <- as.list(1:nclus)
        split.X <- as.list(1:nclus)
        for (k in 1:nclus)
            split.DM[[k]] <- DM[old.clas==k,]            
        # local models computation
        for (k in 1:nclus)
        {   
            nk <- nrow(split.DM[[k]])
            mean.k <- apply(split.DM[[k]],2,mean)# local mean
            sd.k <- sqrt((nk-1)/nk) * apply(split.DM[[k]],2,sd)# local std.dev
            # spliting data matrix for each class
            split.X[[k]] <- scale(split.DM[[k]], center=mean.k, scale=sd.k)
            # calculating outer weights for each class
            out.ws  <- .pls.weights(split.X[[k]], IDM, blocks, modes, scheme, tol, iter)
            w.locals[[k]] <- out.ws[[2]]
            # calculating LV scores for each class
            Y.k <- split.X[[k]] %*% out.ws[[2]]
            # calculating path coefficients for each class
            pathmod <- .pls.paths(IDM, Y.k, plsr)
            path.locals[[k]] <- pathmod[[2]]
            R2.locals[[k]] <- pathmod[[3]][endo==1]
            # calculating loadings and communalities for each class
            loadcomu <- .pls.loads(split.X[[k]], Y.k, blocks)   
            loads.locals[[k]] <- loadcomu[[1]]
            comu.locals[[k]] <- loadcomu[[2]]
            # latent variables for each unit in each local model
            X.k <- scale(DM, center=mean.k, scale=sd.k)
            Y.locals[[k]] <- X.k %*% w.locals[[k]]
            # computation of communality residuals
            out.res <- DM
            for (j in 1:lvs)
            {
                q <- which(blocklist==j) 
                X.hat <- Y.locals[[k]][,j] %*% t(loads.locals[[k]][q])
                out.res[,q] <- (X.k[,q] - X.hat)^2# outer residuals
            }
            outres.locals[[k]] <- out.res
            # computation of inner residuals
            if (sum(endo)!=1)
                Y.hat <- Y.locals[[k]] %*% t(path.locals[[k]][endo==1,])   
            if (sum(endo)==1)
                Y.hat <- Y.locals[[k]] %*% path.locals[[k]][endo==1,]        
            innres.locals[[k]] <- (Y.locals[[k]][,endo==1] - Y.hat)^2
            # computation super normalized residual of outer model
            res.num1 <- outres.locals[[k]] %*% diag(1/comu.locals[[k]],mvs,mvs)
            supres.outer <- rowSums(res.num1) / (sum(rowSums(res.num1))/(N-2))
            # computation of super normalized residual of inner model
            res.num2 <- innres.locals[[k]] %*% diag(1/R2.locals[[k]],n.end,n.end)
            supres.inner <- rowSums(res.num2) / (sum(rowSums(res.num2))/(N-2)) 
            # computation of the CM distance
            CM.locals[,k] <- sqrt(supres.outer * supres.inner)
        }
        # allocating the units to their closest class
        new.clas <- old.clas
        for (i in 1:N)
            new.clas[i] <- which(CM.locals[i,]==min(CM.locals[i,]))[1]
        # checking convergence
        dif.clas <- new.clas - old.clas 
        unit.change <- length(which(dif.clas!=0))
        iter.ch <- iter.ch + 1
        # rate of unit change
        if (unit.change/N < stop.crit || iter.ch==iter.max)
            break
        old.clas <- new.clas
        if (any(table(new.clas)<=5)) 
            stop("Too few units: a class with less than 6 units was detected") 
    }

    # ==================== computation of final local models =======================
    DM.cl.rebus <- as.list(1:nclus)    # data matrices for each class
    for (k in 1:nclus)
        DM.cl.rebus[[k]] <- DM[new.clas==k,]
    DM.rebus <- cbind(DM, rebus.class=new.clas)# data matrix with units membership
    redun.locals <- as.list(1:nclus)
    # table of path coefficients for local models
    path.labs <- NULL
    for (j in 1:lvs)
        for (i in j:lvs)
             if (IDM[i,j]==1) 
                 path.labs <- c(path.labs, paste(lvs.names[j],"->",lvs.names[i],sep=""))
    reb.labs <- paste(rep("Class",nclus),1:nclus,sep=".")
    path.rebus <- matrix(NA,length(path.labs),nclus)
    loads.rebus <- matrix(NA,mvs,nclus)
    qual.rebus <- matrix(NA,sum(lvs+2*n.end+1),nclus)
    for (k in 1:nclus)
    {               
        nk <- nrow(DM.cl.rebus[[k]])
        mean.k <- apply(DM.cl.rebus[[k]], 2, mean)# local mean
        sd.k <- sqrt((nk-1)/nk) * apply(DM.cl.rebus[[k]], 2, sd)# local std.dev
        DM.k <- scale(DM.cl.rebus[[k]], center=mean.k, scale=sd.k)        
        out.ws <- .pls.weights(DM.k, IDM, blocks, modes, scheme, tol, iter)
        Y.k <- DM.k %*% out.ws[[2]]
        pathmod <- .pls.paths(IDM, Y.k, plsr)
        path.locals[[k]] <- pathmod[[2]]
        R2.locals[[k]] <- pathmod[[3]][endo==1]
        loadcomu <- .pls.loads(DM.cl.rebus[[k]], Y.k, blocks)    
        loads.locals[[k]] <- loadcomu[[1]]
        comu.locals[[k]] <- loadcomu[[2]]
        redun <- rep(0,mvs)
        aux <- 0
        for (j in 1:lvs)
            if (endo[j]==1)
            {
                aux <- aux + 1
                redun[blocklist==j] <- comu.locals[[k]][blocklist==j] * R2.locals[[k]][aux]
            }
        redun.locals[[k]] <- redun        
        path.rebus[,k] <- as.vector(path.locals[[k]][IDM==1])# path coeffs for local models        
        loads.rebus[,k] <- loads.locals[[k]]# loadings for local models  
        # table of quality indexes
        comu.aveg <- rep(NA,lvs) 
        redun.aveg <- rep(NA,n.end)
        R2.aux <- rep(NA,n.end)
        aux <- 0
        for (j in 1:lvs)
        {
            comu.aveg[j] <- mean(comu.locals[[k]][which(blocklist==j)]) 
            if (endo[j]==1) 
            {
                aux <- aux + 1      
                redun.aveg[aux] <- mean(redun.locals[[k]][which(blocklist==j)])
                R2.aux[aux] <- R2.locals[[k]][aux]
            }
        }
        qual.rebus[1:lvs,k] <- comu.aveg
        qual.rebus[(lvs+1):(lvs+n.end),k] <- redun.aveg
        qual.rebus[(lvs+n.end+1):(lvs+2*n.end),k] <- R2.aux
        qual.rebus[(lvs+2*n.end+1),k] <- sqrt(mean(comu.aveg)*mean(R2.aux))
    }
    gqi <- round(.pls.GQI(pls, new.clas, DM), 4)
    dimnames(path.rebus) <- list(path.labs, reb.labs)
    dimnames(loads.rebus) <- list(mvs.names, reb.labs)
    v1 <- paste(rep("Com",lvs),lvs.names,sep=".")
    v2 <- paste(rep("Red",n.end),lvs.names[endo==1],sep=".")
    v3 <- paste(rep("R2",n.end),lvs.names[endo==1],sep=".")
    dimnames(qual.rebus) <- list(c(v1,v2,v3,"GoF"), reb.labs)
    qual.rebus <- round(qual.rebus,4)
    aux <- list(lvs,n.end,unit.change/N,stop.crit,iter.max,iter.ch,gqi)
    res <- list(path.coef=path.rebus, loadings=loads.rebus, quality=qual.rebus,
                segments=new.clas, origdata.clas=DM.rebus, aux=aux)
    class(res) <- "rebus"
    return(res)
}

