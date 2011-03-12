plspm.fit <-
function(x, inner, outer, modes=NULL, scheme="centroid", scaled=TRUE,
                      tol=0.00001, iter=100)
{
    # =========================== ARGUMENTS ====================================
    # x: a numeric matrix or data.frame containing the manifest variables
    # inner.mat: a square boolean matrix indicating the path relations
    # sets: a list of vectors of indices indicating the manifest variables 
    #       for each block
    # modes: a character vector indicating the measurement type for each 
    #        latent variable. "A" reflective, "B" formative
    # scheme: a character string indicating the inner weighting scheme 
    #         to be used: "factor", "centroid", or "path"
    # scaled: a logical value indicating whether scale data is performed
    # tol: tolerance threshold for calculating outer weights (0.00001)
    # iter: an integer indicating the maximum number of iterations (100)
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
    plsr <- FALSE
    if (mode(tol)!="numeric" || length(tol)!=1 || tol<=0 || tol>0.001) {
        warning("Warning: Invalid argument 'tol'. Default 'tol=0.00001' is used.")   
        tol <- 0.00001
    } 
    if (mode(iter)!="numeric" || length(iter)!=1 || iter<100) {
        warning("Warning: Invalid argument 'iter'. Default 'iter=100' is used.")   
        iter <- 100
    } 

    # ========================== INPUTS SETTING ==========================
    IDM <- inner
    dimnames(IDM) <- list(lvs.names, lvs.names)
    lvs <- nrow(IDM)
    lvs.names <- rownames(IDM)
    blocks <- unlist(lapply(outer, length))
    mvs <- sum(blocks)
    names(blocks) <- lvs.names
    blocklist <- outer
    for (k in 1:length(outer))
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
        stop("Algorithm stops") 
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
    # ========== Stage 3: Measurement loadings and communalities =========
    loadcomu <- .pls.loads(X, Y.lvs, blocks)    
    loads <- loadcomu[[1]]
    comu <- loadcomu[[2]]
    # ========================= Measurement model ========================
    outmod <- as.list(1:lvs)
    for (j in 1:lvs)
    {
        aux <- which(blocklist==j)
        outmod[[j]] <- round(cbind(weights=out.weights[aux], std.loads=loads[aux], 
                             communal=comu[aux]), 4)
    }
    names(outmod) <- lvs.names
    # =========================== Basic Results ==========================
    skem <- switch(scheme, "centroid"="centroid", "factor"="factor", "path"="path")
    model <- list(IDM=IDM, blocks=blocks, scheme=skem, modes=modes, scaled=scaled, 
                  obs=nrow(X), tol=tol, iter=iter, n.iter=out.ws[[3]], outer=outer)
    res <- list(outer.mod=outmod, inner.mod=innmod, latents=Z.lvs, scores=Y.lvs,
               out.weights=out.weights, loadings=loads, path.coefs=Path, r.sqr=R2, 
               data=NULL, model=model)
    class(res) <- c("plspm.fit", "plspm")
    return(res)
}

