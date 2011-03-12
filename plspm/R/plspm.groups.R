plspm.groups <-
function(pls, g, Y=NULL, method="bootstrap", reps=NULL)
{
    # =========================== ARGUMENTS ==============================
    # pls: an object of class "plspm"
    # g: a factor with 2 levels indicating the groups to be compared
    # Y: optional data matrix used when pls$data is null
    # method: the method to be used: "bootstrap", or "permutation"
    # reps: number of bootstrap resamples or number of permutations
    # ==================== Checking function arguments ===================
    if (class(pls)!="plspm") 
        stop("argument 'pls' must be an object of class 'plspm'")
    if (!is.factor(g)) stop("argument 'g' must be a factor")
    ng <- nlevels(g)
    if (ng > 2) stop("argument 'g' must contain only 2 levels") 
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
    if (!is.na(pmatch(method, "bootstrap"))) 
        method <- "bootstrap"
    METHODS <- c("bootstrap", "permutation")
    method <- pmatch(method, METHODS)
    if (is.na(method)) {
        warning("Invalid argument 'method'. Default 'method=bootstrap' is used.")   
        method <- "bootstrap"
    }
    if (is.null(reps) | length(reps)>1) reps<-100
    if (!is.numeric(reps) | floor(reps)<=0) reps<-100
    # ========================== INPUTS SETTING ==========================
    IDM <- pls$model$IDM# Inner Design Matrix
    blocks <- pls$model$blocks# cardinality of blocks
    scheme <- pls$model$scheme# inner weighting scheme
    modes <- pls$model$modes# measurement modes
    scaled <- pls$model$scaled# type of scaling
    plsr <- pls$model$plsr# pls-regression
    tol <- pls$model$tol# tolerance criterion
    iter <- pls$model$tol# max num iterations
    outer <- pls$model$outer
    blocklist <- outer
    for (k in 1:length(blocks))
         blocklist[[k]] <- rep(k,blocks[k])
    blocklist <- unlist(blocklist)
    # data matrix DM
    if (!is.null(pls$data)) {
        DM <- pls$data
        dataset <- TRUE
    } else {         
        dataset <- FALSE
        # building data matrix 'DM'
        DM <- matrix(NA, nrow(pls$latents), sum(blocks))
        for (k in 1:nrow(IDM))
            DM[,which(blocklist==k)] <- as.matrix(Y[,outer[[k]]])
        dimnames(DM) <- list(rownames(pls$latents), names(pls$out.weights))
    }
    lvs <- nrow(IDM)
    lvs.names <- rownames(IDM)
    mvs <- sum(blocks)
    blocklist <- as.list(1:lvs)
    for (j in 1:lvs)
         blocklist[[j]] <- rep(j,blocks[j])
    blocklist <- unlist(blocklist)
    # apply the selected scaling
    if (scaled) {
        sd.X <- sqrt((nrow(DM)-1)/nrow(DM)) * apply(DM, 2, sd)
        X <- scale(DM, scale=sd.X)
    } else {
        X <- scale(DM, scale=FALSE)
    }
    # ====================== Global model estimation =====================
    out.ws <- .pls.weights(X, IDM, blocks, modes, scheme, tol, iter)
    if (is.null(out.ws)) stop("The pls algorithm is non convergent") 
    cor.XY <- cor(X, X%*%out.ws[[2]])
    w.sig <- rep(NA,lvs)
    for (k in 1:lvs) 
         w.sig[k] <- ifelse(sum(sign(cor.XY[which(blocklist==k),k]))<=0,-1,1)
    if (scaled) {
        Y.lvs <- round(X %*% out.ws[[2]] %*% diag(w.sig,lvs,lvs), 4)
    } else   
        Y.lvs <- round(DM %*% out.ws[[2]] %*% diag(w.sig,lvs,lvs), 4)
    dimnames(Y.lvs) <- list(rownames(X), lvs.names)
    # Path coefficients
    pathmod <- .pls.paths(IDM, Y.lvs, plsr)
    innmod <- pathmod[[1]]
    Path.global <- pathmod[[2]]
    R2.global <- pathmod[[3]]
    endo = rowSums(IDM)
    endo[endo!=0] <- 1  # vector indicating endogenous LVs
    path.labs <- NULL
    efs.labs <- NULL
    for (j in 1:lvs)
        for (i in j:lvs)
             if (IDM[i,j]==1) 
                 path.labs <- c(path.labs, paste(lvs.names[j],"->",lvs.names[i],sep=""))
    path.global <- as.vector(Path.global[which(IDM==1)])
    names(path.global) <- path.labs
 
    # ====================== Group1 model estimation =====================
    g1.lab <- levels(g)[1]
    group1 <- which(g==levels(g)[1])
    ng1 <- length(group1)
    # apply the selected scaling
    if (scaled) {
        sd.Xg1 <- sqrt((ng1-1)/ng1) * apply(DM[group1,], 2, sd)
        X.g1 <- scale(DM[group1,], scale=sd.Xg1)
    } else {
        X.g1 <- scale(DM[group1,], scale=FALSE)
    }
    wgs.g1 <- .pls.weights(X.g1, IDM, blocks, modes, scheme, tol, iter)
    if (is.null(wgs.g1)) stop("The algorithm is non convergent") 
    cor.XY <- cor(X.g1, X.g1%*%wgs.g1[[2]])
    w.sig <- rep(NA,lvs)
    for (k in 1:lvs) 
         w.sig[k] <- ifelse(sum(sign(cor.XY[which(blocklist==k),k]))<=0,-1,1)
    if (scaled) {
        Y1.lvs <- round(X.g1 %*% wgs.g1[[2]] %*% diag(w.sig,lvs,lvs), 4)
    } else   
        Y1.lvs <- round(DM[group1,] %*% wgs.g1[[2]] %*% diag(w.sig,lvs,lvs), 4)
    dimnames(Y1.lvs) <- list(rownames(X.g1), lvs.names)
    # Path coefficients 
    pathmod.g1 <- .pls.paths(IDM, Y1.lvs, plsr)
    innmod.g1 <- pathmod.g1[[1]]
    Path.g1 <- pathmod.g1[[2]]
    R2.g1 <- pathmod.g1[[3]]    
    path.g1 <- as.vector(Path.g1[which(IDM==1)])
    names(path.g1) <- path.labs

    # ====================== Group2 model estimation =====================
    g2.lab <- levels(g)[2]
    group2 <- which(g==levels(g)[2])
    ng2 <- length(group2)
    # apply the selected scaling
    if (scaled) {
        sd.Xg2 <- sqrt((ng2-1)/ng2) * apply(DM[group2,], 2, sd)
        X.g2 <- scale(DM[group2,], scale=sd.Xg2)
    } else {
        X.g2 <- scale(DM[group2,], scale=FALSE)
    }
    wgs.g2 <- .pls.weights(X.g2, IDM, blocks, modes, scheme, tol, iter)
    if (is.null(wgs.g2)) stop("The algorithm is non convergent") 
    cor.XY <- cor(X[group2,], X.g2%*%wgs.g2[[2]])
    w.sig <- rep(NA,lvs)
    for (k in 1:lvs) 
        w.sig[k] <- ifelse(sum(sign(cor.XY[which(blocklist==k),k]))<=0,-1,1)
    if (scaled) {
        Y2.lvs <- round(X.g2 %*% wgs.g2[[2]] %*% diag(w.sig,lvs,lvs), 4)
    } else   
        Y2.lvs <- round(DM[group2,] %*% wgs.g2[[2]] %*% diag(w.sig,lvs,lvs), 4)
    dimnames(Y2.lvs) <- list(rownames(X.g2), lvs.names)
    # Path coefficients 
    pathmod.g2 <- .pls.paths(IDM, Y2.lvs, plsr)
    innmod.g2 <- pathmod.g2[[1]]
    Path.g2 <- pathmod.g2[[2]]
    R2.g2 <- pathmod.g2[[3]]    
    path.g2 <- as.vector(Path.g2[which(IDM==1)])
    names(path.g2) <- path.labs

    # ====================== Group Comparison =====================
    dif.orig <- abs(path.g1-path.g2)
    nb <- round(reps)

    if (method==1)   # bootstrap
    {
        BG1 <- matrix(0, nb, sum(IDM))
        BG2 <- BG1
        for (i in 1:nb)
        {
            samg1 <- sample(group1, ng1, replace=TRUE) 
            samg2 <- sample(group2, ng2, replace=TRUE)
            # apply the selected scaling
            if (scaled) {
                sd.Xg1 <- sqrt((ng1-1)/ng1) * apply(DM[samg1,], 2, sd)
                sd.Xg2 <- sqrt((ng2-1)/ng2) * apply(DM[samg2,], 2, sd)
                X.g1 <- scale(DM[samg1,], scale=sd.Xg1)
                X.g2 <- scale(DM[samg2,], scale=sd.Xg2)
            } else {
                X.g1 <- scale(DM[samg1,], scale=FALSE)
                X.g2 <- scale(DM[samg2,], scale=FALSE)
            }
            wgs.g1 <- .pls.weights(X.g1, IDM, blocks, modes, scheme, tol, iter)
            wgs.g2 <- .pls.weights(X.g2, IDM, blocks, modes, scheme, tol, iter)
            if (is.null(wgs.g1)) stop("Non convergence in bootstrap samples") 
            if (is.null(wgs.g2)) stop("Non convergence in bootstrap samples") 
            cor.XY <- cor(X.g1, X.g1%*%wgs.g1[[2]])
            w.sig <- rep(NA,lvs)
            for (k in 1:lvs) 
                 w.sig[k] <- ifelse(sum(sign(cor.XY[which(blocklist==k),k]))<=0,-1,1)
            if (scaled) {
                Y1.lvs <- X.g1 %*% wgs.g1[[2]] %*% diag(w.sig,lvs,lvs)
            } else   
                Y1.lvs <- DM[samg1,] %*% wgs.g1[[2]] %*% diag(w.sig,lvs,lvs)
            cor.XY <- cor(X.g2, X.g2%*%wgs.g2[[2]])
            w.sig <- rep(NA,lvs)
            for (k in 1:lvs) 
                 w.sig[k] <- ifelse(sum(sign(cor.XY[which(blocklist==k),k]))<=0,-1,1)
            if (scaled) { 
                Y2.lvs <- X.g2 %*% wgs.g2[[2]] %*% diag(w.sig,lvs,lvs)
            } else   
                Y2.lvs <- DM[samg2,] %*% wgs.g2[[2]] %*% diag(w.sig,lvs,lvs)
            pathmod.g1 <- .pls.paths(IDM, Y1.lvs, plsr)
            paths.g1 <- pathmod.g1[[2]]    
            pathmod.g2 <- .pls.paths(IDM, Y2.lvs, plsr)
            paths.g2 <- pathmod.g2[[2]]
            BG1[i,] <- as.vector(paths.g1[which(IDM==1)])
            BG2[i,] <- as.vector(paths.g2[which(IDM==1)])
        }    
        path.difs <- abs(apply(BG1,2,mean) - apply(BG2,2,mean))
        SE1 <- apply(BG1, 2, var)
        SE2 <- apply(BG2, 2, var)
        names(path.global) <- path.labs
        t.stat <- rep(NA, sum(IDM))
        k1 <- ((ng1-1)^2)/(ng1+ng2-2)
        k2 <- ((ng2-1)^2)/(ng1+ng2-2)
        k3 <- sqrt(1/ng1 + 1/ng2)
        for (i in 1:sum(IDM))         
             t.stat[i] <- path.difs[i] / (sqrt(k1*SE1[i]+k2*SE2[i]) * k3)        
        p.val <- pt(t.stat, ng1+ng2-2, lower.tail=FALSE)
        signi.path <- rep("no",length(p.val))
        signi.path[p.val<0.05] <- "yes"
        res.path <- round(cbind(path.global, path.g1, path.g2, dif.orig, 
                        t.stat, df=rep(ng1+ng2-2,sum(IDM)), p.val), 4)
        res <- data.frame(res.path, signi.path)
        colnames(res) <- c("global", paste(rep("group",2),levels(g),sep="."), 
                           "diff.abs", "t.stat", "deg.fr", "p.value", "sig.05") 
    } else
    {
        dif.perm <- matrix(0, nb, sum(IDM))
        for (i in 1:nb)# multigroup permutation
        {
            permu <- sample(1:(ng1+ng2), ng1+ng2)
            samg1 <- permu[1:ng1]
            samg2 <- permu[(ng1+1):(ng1+ng2)]
            # apply the selected scaling
            if (scaled) {
                sd.Xg1 <- sqrt((ng1-1)/ng1) * apply(DM[samg1,], 2, sd)
                sd.Xg2 <- sqrt((ng2-1)/ng2) * apply(DM[samg2,], 2, sd)
                X.g1 <- scale(DM[samg1,], scale=sd.Xg1)
                X.g2 <- scale(DM[samg2,], scale=sd.Xg2)
            } else {
                X.g1 <- scale(DM[samg1,], scale=FALSE)
                X.g2 <- scale(DM[samg2,], scale=FALSE)
            }
            wgs.g1 <- .pls.weights(X.g1, IDM, blocks, modes, scheme, tol, iter)
            wgs.g2 <- .pls.weights(X.g2, IDM, blocks, modes, scheme, tol, iter)
            if (is.null(wgs.g1)) stop("Non convergence in bootstrap samples") 
            if (is.null(wgs.g2)) stop("Non convergence in bootstrap samples") 
            cor.XY <- cor(X.g1, X.g1%*%wgs.g1[[2]])
            w.sig <- rep(NA,lvs)
            for (k in 1:lvs) 
                 w.sig[k] <- ifelse(sum(sign(cor.XY[which(blocklist==k),k]))<=0,-1,1)
            if (scaled) {
                Y1.lvs <- X.g1 %*% wgs.g1[[2]] %*% diag(w.sig,lvs,lvs)
            } else   
                Y1.lvs <- DM[samg1,] %*% wgs.g1[[2]] %*% diag(w.sig,lvs,lvs)
            cor.XY <- cor(X.g2, X.g2%*%wgs.g2[[2]])
            w.sig <- rep(NA,lvs)
            for (k in 1:lvs) 
                 w.sig[k] <- ifelse(sum(sign(cor.XY[which(blocklist==k),k]))<=0,-1,1)
            if (scaled) {
                Y2.lvs <- X.g2 %*% wgs.g2[[2]] %*% diag(w.sig,lvs,lvs)
            } else   
                Y2.lvs <- DM[samg2,] %*% wgs.g2[[2]] %*% diag(w.sig,lvs,lvs)
            pathmod.g1 <- .pls.paths(IDM, Y1.lvs, plsr)
            paths.g1 <- pathmod.g1[[2]]    
            pathmod.g2 <- .pls.paths(IDM, Y2.lvs, plsr)
            paths.g2 <- pathmod.g2[[2]]
            pp1 <- as.vector(paths.g1[which(IDM==1)])
            pp2 <- as.vector(paths.g2[which(IDM==1)])
            dif.perm[i,] <- abs(pp1 - pp2)
        }   
        s.perm <- dif.orig 
        for (i in 1:sum(IDM))         
            s.perm[i] <- length(which(dif.orig[i]<dif.perm[,i])) + 1
        p.val <- (1/(nb+1))*s.perm 
        signi.path <- rep("no",length(p.val))
        signi.path[p.val<0.05] <- "yes"
        res.path <- round(cbind(path.global, path.g1, path.g2, dif.orig, p.val), 4)
        res <- data.frame(res.path, signi.path)
        colnames(res) <- c("global", paste(rep("group",2),levels(g),sep="."), 
                             "diff.abs", "p.value", "sig.05") 
    }
    met <- switch(method, "1"="bootstrap", "2"="permutation")
    settings <- c(scaled=scaled, scheme=scheme, method=met)
    res <- list(test=res, global=innmod, group1=innmod.g1, group2=innmod.g2, 
                settings=settings, reps=reps)
    class(res) <- "plspm.groups"
    return(res)
}

