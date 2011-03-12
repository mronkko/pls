res.clus <-
function(pls, Y=NULL)
{
    # ========================= res.clus function ========================
    # Function to calculate communality and structural residuals that
    # will be used for the Response-Based Units Segmentation (REBUS)
    # in the "rebus" function
    # =========================== ARGUMENTS ==============================
    # pls: object of class "plspm"
    # Y: optional data matrix used when pls$data is null
    ## pls$model <- list(IDM, blocks, scheme, modes, scaled, boot.val, 
    ##                    plsr, obs, br, tol, iter, n.iter, outer)

    # ==================== Checking function arguments ===================
    if (class(pls)!="plspm") 
        stop("An object of class 'plspm' was expected")
    if (any(pls$model$modes!="A"))# checking reflective modes
        stop("REBUS only works for reflective modes")
    if (!pls$model$scaled)# checking scaled data
        stop("REBUS only works with scaled='TRUE'")
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
    scheme <- pls$model$scheme# inner weighting scheme
    modes <- pls$model$modes# measurement modes
    scaled <- pls$model$scaled# type of scaling
    plsr <- pls$model$plsr# pls-regression
    tol <- pls$model$tol# tolerance criterion
    iter <- pls$model$iter# max num iterations
    outer <- pls$model$outer
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
    lvs <- nrow(IDM)
    lvs.names <- rownames(IDM)
    mvs <- sum(blocks)
    # data scaling (standardized data)
    sd.X <- sqrt((nrow(DM)-1)/nrow(DM)) * apply(DM, 2, sd)
    X <- scale(DM, scale=sd.X)
   
    # ====================== computation of residuals =====================
    Y.lvs <- pls$latents# recovering LV scores from pls
    loads <- pls$loadings# recovering loadings from pls
    PaCo <- pls$path.coefs# recovering path coeffs from pls
    endo <- rowSums(IDM)
    endo[endo!=0] <- 1  # indicator of endogenous LVs
    out.res <- DM# matrix for storing outer resids
    inn.res <- Y.lvs[,endo==1]# matrix for storing inner resids
    # computation of outer residuals
    for (j in 1:lvs)
    {
        q <- which(blocklist==j) 
        X.hat <- Y.lvs[,j] %*% t(loads[q])
        out.res[,q] <- X[,q] - X.hat# outer residuals
    }
    # computationi of inner residuals
    if (sum(endo)!=1)# more than 1 endogenous LV
        Y.hat <- Y.lvs %*% t(PaCo[endo==1,])        
    if (sum(endo)==1)# only 1 endogenous LV
        Y.hat <- Y.lvs %*% PaCo[endo==1,]        
    inn.res <- Y.lvs[,endo==1] - Y.hat# inner residuals
    
    # ====================== cluster analysis =====================
    # hierarchical cluster analysis with Ward method using function "hcluster"
    res <- cbind(out.res, inn.res)    
    res.clus <- hcluster(res, method="euclidean", diag=FALSE, upper=FALSE,
     link="ward", members=NULL, nbproc=2, doubleprecision=TRUE)
    # plot of the dendrogram
    plot(res.clus, main=c("REBUS", "Cluster Dendrogram of Outer and Inner Residuals"),
         hang=-1, cex.main=.9, cex.axis=.5, xlab="Hierarchical Clustering", 
         sub="Ward method", labels=FALSE)
    return(res.clus)
}

