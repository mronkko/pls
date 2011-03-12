local.models <-
function(pls, y, Y=NULL)
{
    # ======================== local.models function ======================
    # Function to calculate PLS-PM for global and local models
    # using a "rebus" object or a vector (with memberships or categories)
    # =========================== arguments ===============================
    # pls: object of class "plspm"
    # y: can be an object of class "rebus", a numeric vector, or a factor
    # Y: optional data matrix or dataframe

    # ==================== Checking function arguments ====================
    if (class(pls)!="plspm") 
        stop("argument 'pls' must be an object of class 'plspm'")
    if (!is.element(class(y), c("rebus","integer","factor")))   
        stop("argument 'y' must be of class 'rebus', 'integer' or 'factor'")
    if (class(y)=="rebus") {
        if (length(y$segments)!=nrow(pls$latents))
            stop("arguments 'pls' and 'y' are incompatible")
    } else {
        if (length(y)!=nrow(pls$latents))
            stop("arguments 'pls' and 'y' are incompatible")
    }
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

    # ========================== INPUTS SETTING ==========================
    IDM <- pls$model$IDM# Inner Design Matrix
    blocks <- pls$model$blocks# cardinality of blocks
    modes <- pls$model$modes# measurement modes    
    plsr <- FALSE 
    tol <- pls$model$tol
    iter <- pls$model$iter
    scheme <- pls$model$scheme
    scaled <- pls$model$scaled
    tol <- pls$model$tol
    iter <- pls$model$iter
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
    endo <- rowSums(IDM)
    endo[endo!=0] <- 1  
    end.ind <- cumsum(blocks)
    ini.ind <- cumsum(blocks) - blocks + 1
    new.sets <- as.list(1:lvs)
    for (j in 1:lvs)
        new.sets[[j]] <- ini.ind[j]:end.ind[j]
    if (class(y)=="rebus") {
        segments <- as.factor(y$segments)
    } else {
        segments <- as.factor(y)
    }
    n.clus <- length(table(segments))

    # ============ final models computation (global and local models) ============
    skem <- switch(scheme, "centroid"="centroid", "factor"="factor")
    final.mod <- as.list(1:(n.clus+1))# final plspm models
    for (k in 1:(n.clus+1))
    {
        if (k==1) {
            # global model
            X <- DM
            final.mod[[1]] <- plspm(X, IDM, new.sets, modes, skem, scaled, 
                              tol=tol, iter=iter, dataset=dataset)
        } else
        {
            units.k <- which(segments==levels(segments)[k-1])
            # local models
            X.k <- DM[units.k,]
            final.mod[[k]] <- plspm(X.k, IDM, new.sets, modes, skem, scaled, 
                              tol=tol, iter=iter, dataset=dataset)
        }
    }
    names(final.mod) <- c("glob.model",paste(rep("loc.model",n.clus), 1:n.clus, sep="."))
    class(final.mod) <- "local.models"
    return(final.mod)
}

