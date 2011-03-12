plspm <-
function(x, inner, outer, modes=NULL, scheme="centroid", scaled=TRUE,
            boot.val=FALSE, br=NULL, plsr=FALSE, tol=0.00001, iter=100, dataset=TRUE)
{
    # =========================== ARGUMENTS ====================================
    # x: a numeric matrix or data.frame containing the manifest variables
    # inner: a square boolean matrix representing the inner model
    # outer: a list representing the outer model (list of vectors)
    #        from matrix x with manifest variables for each block
    # modes: a character vector indicating the measurement type for each 
    #        latent variable. "A" reflective, "B" formative
    # scheme: a character string indicating the inner weighting scheme 
    #         to be used: "factor", "centroid", or "path"
    # scaled: a logical value indicating whether scale data is performed
    # boot.val:a logical value indicating whether bootstrap validation is done 
    # br: an integer indicating the number of bootstraps resamples, used 
    #     only when boot.val=TRUE, (100 <= br <= 1000)
    # plsr: a logical value for calculating path coeffs by pls-regression
    # tol: tolerance threshold for calculating outer weights (0.00001)
    # iter: an integer indicating the maximum number of iterations (100)
    # dataset: a logical value indicating whether the data should be retrieved
    # ===========================================================================

    # ==================== Checking function arguments ===================
    if (!is.matrix(x) && !is.data.frame(x))
        stop("Invalid object 'x'. Must be a numeric matrix or data frame.")
    if (is.null(rownames(x))) 
        rownames(x) <- 1:nrow(x)
    if (is.null(colnames(x))) 
        colnames(x) <- paste("MV", 1:ncol(x), sep="")
    if (!is.matrix(inner))
        stop("Invalid argument 'inner'. Must be a matrix.")
    if (nrow(inner)!=ncol(inner))
        stop("Invalid argument 'inner'. Must be a square matrix.")
    for (j in 1:ncol(inner))
        for (i in 1:nrow(inner)) {
            if (i<=j)
                if (inner[i,j]!=0) 
                    stop("argument 'inner' must be a lower triangular matrix")
            if (length(intersect(inner[i,j], c(1,0)))==0)
                stop("elements in 'inner' must be '1' or '0'")
        }
    if (is.null(dimnames(inner)))
        lvs.names <- paste("LV", 1:ncol(inner), sep="")
    if (!is.null(rownames(inner)))
        lvs.names <- rownames(inner)
    if (!is.null(colnames(inner)))
        lvs.names <- colnames(inner)
    if (!is.list(outer))
        stop("Invalid argument 'outer'. Must be a list.")
    if (length(outer) != nrow(inner))
        stop("Number of rows of 'inner' does not coincide with length of 'outer'.")
    if (is.null(modes)) {
        modes <- rep("A",length(outer))
        warning("Argument 'modes' missing. Default reflective 'modes' is used.")
    }
    if (length(outer) != length(modes)) {
        warning("Warning: Invalid length of 'modes'. Default reflective 'modes' is used.")
        modes <- rep("A", length(outer))
    }
    for (i in 1:length(modes))
        if (modes[i]!="A" && modes[i]!="B") modes[i]<-"A"
    if (!is.na(pmatch(scheme, "centroid"))) 
        scheme <- "centroid"
    SCHEMES <- c("centroid", "factor", "path")
    scheme <- pmatch(scheme, SCHEMES)
    if (is.na(scheme)) {
        warning("Warning: Invalid argument 'scheme'. Default 'scheme=centroid' is used.")   
        scheme <- "centroid"
    }
    if (!is.logical(scaled)) {
        warning("Warning: Invalid argument 'scaled'. Default 'scaled=TRUE' is used.")
        scaled <- TRUE
    }
    if (!is.logical(boot.val)) {
        warning("Warning: Invalid argument 'boot.val'. No bootstrap validation is done.")
        boot.val <- FALSE
    }   
    if (boot.val) {
        if (!is.null(br)) {        
            if (mode(br)!="numeric" || length(br)!=1 || (br%%1)!=0 ||
                br<100 || br>1000) {
                warning("Warning: Invalid argument 'br'. Default 'br=100' is used.")   
                br <- 100
            } 
        } else
            br <- 100
    }
    if (!is.logical(plsr)) plsr<-FALSE
    if (mode(tol)!="numeric" || length(tol)!=1 || tol<=0 || tol>0.001) {
        warning("Warning: Invalid argument 'tol'. Default 'tol=0.00001' is used.")   
        tol <- 0.00001
    } 
    if (mode(iter)!="numeric" || length(iter)!=1 || iter<100) {
        warning("Warning: Invalid argument 'iter'. Default 'iter=100' is used.")   
        iter <- 100
    } 
    if (!is.logical(dataset)) 
        dataset <- TRUE

    # ========================== INPUTS SETTING ==========================
    IDM <- inner
    dimnames(IDM) <- list(lvs.names, lvs.names)
    lvs <- nrow(IDM)
    lvs.names <- rownames(IDM)
    blocks <- unlist(lapply(outer, length))
    mvs <- sum(blocks)
    names(blocks) <- lvs.names
    blocklist <- outer
    for (k in 1:length(blocks))
         blocklist[[k]] <- rep(k,blocks[k])
    blocklist <- unlist(blocklist)
    Mode <- modes
    Mode[modes=="A"] <- "Reflective"
    Mode[modes=="B"] <- "Formative"   
    # building data matrix 'DM'
    DM <- matrix(NA, nrow(x), mvs)
    mvs.names <- rep(NA, mvs)
    for (k in 1:lvs)
    {        
        DM[,which(blocklist==k)] <- as.matrix(x[,outer[[k]]])
        mvs.names[which(blocklist==k)] <- colnames(x)[outer[[k]]]
    }
    dimnames(DM) <- list(rownames(x), mvs.names)
    # apply the selected scaling
    if (scaled) {
        sd.X <- sqrt((nrow(DM)-1)/nrow(DM)) * apply(DM, 2, sd)
        X <- scale(DM, scale=sd.X)
    } else {
        X <- scale(DM, scale=FALSE)
    }
    dimnames(X) <- list(rownames(x), mvs.names)

    # ==================== Stage 1: Iterative procedure ==================
    out.ws <- .pls.weights(X, IDM, blocks, modes, scheme, tol, iter)
    if (is.null(out.ws)) {
        print(paste("Iterative process is non-convergent with 'iter'=", iter, " and 'tol'=", tol, sep=""))
        message("Algorithm stops") 
        stop("")
    }
    out.weights <- out.ws[[1]]
    cor.XY <- cor(X, X%*%out.ws[[2]])
    w.sig <- rep(NA,lvs)
    for (k in 1:lvs) 
         w.sig[k] <- ifelse(sum(sign(cor.XY[which(blocklist==k),k]))<=0,-1,1)
    Z.lvs <- X %*% out.ws[[2]] %*% diag(w.sig,lvs,lvs)
    Y.lvs <- Z.lvs
    if (!scaled) 
        Y.lvs <- DM %*% out.ws[[2]] %*% diag(w.sig,lvs,lvs)
    dimnames(Y.lvs) <- list(rownames(X), lvs.names)
    dimnames(Z.lvs) <- list(rownames(X), lvs.names)
    # ============ Stage 2: Path coefficients and total effects ==========
    pathmod <- .pls.paths(IDM, Y.lvs, plsr)
    innmod <- pathmod[[1]]
    Path <- pathmod[[2]]
    R2 <- pathmod[[3]]
    Path.efs <- .pls.efects(Path)
    # ========== Stage 3: Measurement loadings and communalities =========
    loadcomu <- .pls.loads(X, Y.lvs, blocks)    
    loads <- loadcomu[[1]]
    comu <- loadcomu[[2]]
    endo <- rowSums(IDM)
    endo[endo!=0] <- 1  
    redun <- rep(0,mvs)
    for (j in 1:lvs)
        if (endo[j]==1)
            redun[blocklist==j] <- comu[blocklist==j] * R2[j]
    # ========================= Measurement model ========================
    outcor <- as.list(1:lvs)
    outmod <- as.list(1:lvs)
    for (j in 1:lvs)
    {
        aux <- which(blocklist==j)
        outmod[[j]] <- round(cbind(weights=out.weights[aux], std.loads=loads[aux], 
                             communal=comu[aux], redundan=redun[aux]), 4)
        outcor[[j]] <- round(cor(DM[,aux], Y.lvs), 4)
    }
    names(outmod) <- lvs.names
    names(outcor) <- lvs.names  
    # ======================== Unidimensionality =========================
    unidim <- .pls.unidim(DM, blocks, modes)
    # ======================== Summary Inner model =======================
    exo.endo <- rowSums(IDM)
    exo.endo[rowSums(IDM)==0] <- "Exogen"
    exo.endo[rowSums(IDM)!=0] <- "Endogen"
    av.comu <- rep(0,lvs)   # average communality
    av.redu <- rep(0,lvs)   # average redundancy
    ave <- rep(0, lvs)      # average variance extracted
    for (k in 1:lvs)
    {
        av.comu[k] <- mean(comu[which(blocklist==k)])
        av.redu[k] <- mean(redun[which(blocklist==k)])
        if (modes[k]=="A")
        {
            ave.num <- sum(comu[which(blocklist==k)])
            ave.denom <- sum(comu[which(blocklist==k)]) + sum(1-(comu[which(blocklist==k)]))
            ave[k] <- round(ave.num / ave.denom, 3)
        }
    }
    names(ave) <- lvs.names
    innsum = data.frame(LV.Type=exo.endo, Measure=abbreviate(Mode,5), MVs=blocks, 
           R.square=round(R2,4), Av.Commu=round(av.comu,4), Av.Redun=round(av.redu,4), AVE=ave)
    rownames(innsum) <- lvs.names
    # ============================ GoF Indexes ===========================
    gof <- .pls.GOF(DM, IDM, blocks, comu, unidim, R2)
    # ============================= Results ==============================
    skem <- switch(scheme, "centroid"="centroid", "factor"="factor", "path"="path")
    model <- list(IDM=IDM, blocks=blocks, scheme=skem, modes=modes, scaled=scaled, 
                  boot.val=boot.val, plsr=plsr, obs=nrow(X), br=br, 
                  tol=tol, iter=iter, n.iter=out.ws[[3]], outer=outer)
    if (dataset) data=DM else data=NULL
    res <- list(outer.mod=outmod, inner.mod=innmod, latents=Z.lvs, scores=Y.lvs,
               out.weights=out.weights, loadings=loads, path.coefs=Path, r.sqr=R2,
               outer.cor=outcor, inner.sum=innsum, effects=Path.efs, unidim=unidim, gof=gof, 
               data=data, model=model)
    # ========================= Bootstrap Validation =========================
    if (boot.val) 
    {
        if (nrow(X) <= 10) {
            warning("Bootstrapping stopped: very few cases.") 
        } else 
        { 
            n.efs <- nrow(Path.efs)
            res.boot <- .pls.boot(DM, IDM, blocks, modes, scheme, scaled, br, plsr, tol, iter)
        }
        res <- list(outer.mod=outmod, inner.mod=innmod, latents=Z.lvs, scores=Y.lvs,
                   out.weights=out.weights, loadings=loads, path.coefs=Path, r.sqr=R2,
                   outer.cor=outcor, inner.sum=innsum, effects=Path.efs, unidim=unidim, gof=gof, 
                   boot=res.boot, data=data, model=model)
    } # end 'if' bootstrapping
    # --------------------------------------------------------------------
    class(res) <- "plspm"
    return(res)
}

