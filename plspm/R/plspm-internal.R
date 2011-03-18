.cutree.order <-
function(hclu, k=NULL, h=NULL)
{  
    coupe <- cutree(hclu, k=k, h=h)
    coupe.or <- coupe[hclu$order]
    coupe.out<- rep(NA,length(coupe))
    j <- 1 
    k <- coupe.or[1]
    for (i in 1:length(coupe))
    {
        if (coupe.or[i]==k) next
        else {
            coupe.out[which(coupe==k)] <- j
            j <- j + 1
            k <- coupe.or[i]
        }
    }
    coupe.out[is.na(coupe.out)] <- j
    names(coupe.out) <- names(coupe)
    coupe.out
}

.innerplot <-
function(path.coefs, arr.pos=arr.pos, box.prop=box.prop, 
                       box.cex=box.cex, cex.txt=cex.txt)
{
    # plot of inner model
    MPC <- path.coefs    # matrix of path coefficients
    AM.col <- MPC    # arrow matrix colors
    AM.col[MPC<0] <- "red"# negative path coeffs in red
    AM.col[MPC>0] <- "blue"   # positive path coeffs in blue
    names <- rownames(MPC)# vector of names
    dev.new()
    plotmat(round(MPC,4), pos=NULL, curve=0, name=names, lwd=1, box.lwd=1, cex.txt=cex.txt,
        box.type="circle", box.prop=box.prop, box.cex=box.cex, arr.type="triangle", arr.pos=arr.pos,
        shadow.size=0.01, prefix="", arr.lcol=AM.col, arr.col=AM.col, arr.width=.2,
        main=c("Inner Model","Path Coefficients"))
}

.loadingsplot <-
function(IDM, modes, blocks, loadings, arr.pos=arr.pos,
                     box.prop=box.prop, box.cex=box.cex, cex.txt=cex.txt, newdev=newdev)
{
    lvs <- nrow(IDM)
    ini.vec <- cumsum(blocks) - blocks + 1
    end.vec <- cumsum(blocks) 
    ## plot of loadings
    for (k in 1:lvs)
    {
        num.mvs <- blocks[[k]]
        names.mvs <- names(loadings)[ini.vec[k]:end.vec[k]]
        names.mvs <- c(names.mvs,rownames(IDM)[k])
        box.types <- c(rep("rect",num.mvs),"ellipse")
        ML <- matrix(0,num.mvs+1,num.mvs+1)
        ML.col <- ML
        ML[num.mvs+1,] <- c(loadings[ini.vec[k]:end.vec[k]],0)
        ML.col[ML<0] <- "red"# negative loadinsg in red
        ML.col[ML>0] <- "blue"   # positive loadings in blue
        if (newdev)
            dev.new()
        if (modes[k]=="A") {   # mode "A"
            plotmat(round(t(ML),4), curve=0, name=names.mvs, lwd=1, box.type=box.types, arr.width=.1,
               arr.pos=arr.pos, arr.lcol=t(ML.col), arr.col=t(ML.col), box.prop=box.prop, box.cex=box.cex,
               cex.txt=cex.txt, main=c(paste(rownames(IDM)[k]),"loadings"))
        } else {   # mode "B"
                plotmat(round(ML,4), curve=0, name=names.mvs, lwd=1, box.type=box.types, arr.width=.1,
                   arr.pos=arr.pos, arr.lcol=ML.col, arr.col=ML.col, box.prop=box.prop, box.cex=box.cex,
                   cex.txt=cex.txt, main=c(paste(rownames(IDM)[k]),"loadings"))
        }
    }
}

.path.scheme <-
function(IDM, Y)
{
    lvs <- nrow(IDM)
    E <- IDM
    for (k in 1:lvs) 
    {
        if (length(which(IDM[k,]==1)) > 0)
            E[which(IDM[k,]==1),k] <- lm(Y[,k]~Y[,which(IDM[k,]==1)]-1)$coef
        if (length(which(IDM[,k]==1)) > 0)
            E[which(IDM[,k]==1),k] <- cor(Y[,k], Y[,which(IDM[,k]==1)])
    }                 
    return(E)
}

#
# DM: Data matrix
# IDM: Inner model
# blocks: outer model
# modes: a character vector indicating the measurement type for each 
#        latent variable. "A" reflective, "B" formative
# scheme: a character string indicating the inner weighting scheme 
#         to be used: "factor", "centroid", or "path"
# scaled: a logical value indicating whether scale data is performed
# br: an integer indicating the number of bootstraps resamples, used 
# plsr: a logical value for calculating path coeffs by pls-regression
# tol: tolerance threshold for calculating outer weights (0.00001)
# iter: an integer indicating the maximum number of iterations (100)
# orig.weights: original weights. Used for sign corrections

.pls.boot <-
function(DM, IDM, blocks, modes, scheme, scaled, br, plsr, tol, iter,orig.weights)
{

	signChangeCorrections<-c("Standard","IndividualSignChanges","ConstructLevelChanges")
	
    lvs <- nrow(IDM)
    lvs.names <- rownames(IDM)
    mvs <- ncol(DM)
    mvs.names <- colnames(DM)
    blocklist <- as.list(1:lvs)
    for (j in 1:lvs)
         blocklist[[j]] <- rep(j,blocks[j])
    blocklist <- unlist(blocklist)
    endo <- rowSums(IDM)
    endo[endo!=0] <- 1    
    bootnum <- br
    # scaling data
    if (scaled) {
        sd.X <- sqrt((nrow(DM)-1)/nrow(DM)) * apply(DM, 2, sd)
        X <- scale(DM, scale=sd.X)
    } else {
        X <- scale(DM, scale=FALSE)
    }
    colnames(X) <- mvs.names
    # =============== computation of the original plspm model ================
    out.ws <- .pls.weights(X, IDM, blocks, modes, scheme, tol, iter)
    wgs.orig <- out.ws[[1]]
    cor.XY <- cor(X, X%*%out.ws[[2]])
    w.sig <- rep(NA,lvs)
    for (k in 1:lvs) 
         w.sig[k] <- ifelse(sum(sign(cor.XY[which(blocklist==k),k]))<=0,-1,1)
    Y.lvs <- X %*% out.ws[[2]] %*% diag(w.sig,lvs,lvs)
    pathmod <- .pls.paths(IDM, Y.lvs, plsr)
    Path <- pathmod[[2]]
    path.orig <- as.vector(Path[which(IDM==1)])
    r2.orig <- pathmod[[3]][which(endo==1)]
    Path.efs <- .pls.efects(Path)
    loadcomu <- .pls.loads(X, Y.lvs, blocks)    
    loads.orig <- loadcomu[[1]]
    # ========================= Bootstrap Validation =========================
    path.labs <- NULL
    efs.labs <- NULL
    for (j in 1:lvs)
        for (i in j:lvs)
             if (IDM[i,j]==1) 
                 path.labs <- c(path.labs, paste(lvs.names[j],"->",lvs.names[i],sep=""))
    
    
    WEIGS <- matrix(NA, bootnum, mvs)
    LOADS <- matrix(NA, bootnum, mvs)
    PATHS <- matrix(NA, bootnum, sum(IDM))
    TOEFS <- matrix(NA, bootnum, nrow(Path.efs))
    RSQRS <- matrix(NA, bootnum, sum(endo))
	
	boot.all<-list()
	for(i in 1:length(signChangeCorrections)){ 
		boot.all[[signChangeCorrections[i]]]<-list(WEIGS=WEIGS,LOADS=LOADS,PATHS=PATHS,TOEFS=TOEFS,RSQRS=RSQRS)
	}
	#
	# The main bootsrapping loop
	#
    i <- 1
    W <- NULL
    while (i <= bootnum)
    {
        boot.obs <- sample.int(nrow(X), size=nrow(X), replace=TRUE)
        DM.boot <- DM[boot.obs,]
        # scaling boot sample
        if (scaled) {
            sd.XB <- sqrt((nrow(DM.boot)-1)/nrow(DM.boot)) * apply(DM.boot, 2, sd)
            X.boot <- scale(DM.boot, scale=sd.XB)
        } else {
            X.boot <- scale(DM.boot, scale=FALSE)
        }
        colnames(X.boot) <- mvs.names

        # calculating boot model parameters. Use the weights from previous iteration as starting
        # values. If this is the first iteration, set startign values to null

        w.boot <- .pls.weights(X.boot, IDM, blocks, modes, scheme, tol, iter,W)
        if (is.null(w.boot)) {
            i <- i - 1
            next
        }
        
        #
        # Apply the three different sign change corrections and store the results
        #

        W <- w.boot[[2]]

		for(i in 1:length(signChangeCorrections)){ 
			
			currentResults<-boot.all[[signChangeCorrections[i]]]
    		
    		corrected.w<-w.boot[[1]]
    		corrected.W<-w.boot[[2]]
    		if(signChangeCorrections[i]=="IndividualSignChanges"){
    			correctionVector=sign(corrected.w*orig.weights)
    			corrected.w<-correctionVector*corrected.w
    			corrected.W<-correctionVector*corrected.W
    		}
    		else if(signChangeCorrections[i]=="ConstructLevelChanges"){
    			bitmatrix<-corrected.W!=0
    			flipConstructs<-colSums(bitmatrix*orig.weights-corrected.W)>colSums(bitmatrix*orig.weights+corrected.W)
				flipMatrix<-matrix(flipConstructs,nrow = nrow(corrected.W), ncol = ncol(corrected.W), byrow =TRUE)
    			corrected.W[flipMatrix]<- - corrected.W[flipMatrix]
    			corrected.w<-rowSums(corrected.W)
    		}
    		
	        currentResults[["WEIGS"]][i,] <- corrected.w
	        Y.boot <- X.boot %*% corrected.W
	        pathmod <- .pls.paths(IDM, Y.boot, plsr)
	        P.boot <- pathmod[[2]]
	        Toef.boot <- .pls.efects(P.boot)
	        currentResults[["PATHS"]][i,] <- as.vector(P.boot[which(IDM==1)])
	        currentResults[["TOEFS"]][i,] <- Toef.boot[,4]
	        currentResults[["RSQRS"]][i,] <- pathmod[[3]][which(endo==1)]
	        l.boot <- .pls.loads(X.boot, Y.boot, blocks)    
	        currentResults[["LOADS"]][i,] <- l.boot[[1]]
	        
	        boot.all[[signChangeCorrections[i]]]<-currentResults
	    }
        i <- i + 1
    }
    
    #
    # End of the main bootsrapping loop
    #
    
    # Outer weights

	res.boot<-list()
	
	for(i in 1:length(signChangeCorrections)){ 
		
		currentResults<-boot.all[[signChangeCorrections[i]]]
		WEIGS<-currentResults$WEIGS
		LOADS<-currentResults$LOADS
		PATHS<-currentResults$PATHS
		TOEFS<-currentResults$TOEFS
		RSQRS<-currentResults$RSQRS

	    colnames(WEIGS) <- mvs.names
	    WB <- data.frame(Original = wgs.orig, Mean.Boot = apply(WEIGS, 2, mean), 
	        Std.Error = apply(WEIGS, 2, sd), perc.05 = apply(WEIGS, 2, function(x) quantile(x, 0.05)),
	        perc.95 = apply(WEIGS, 2, function(x) quantile(x, 0.95)))
	    colnames(LOADS) <- mvs.names
	    LB <- data.frame(Original = loads.orig, Mean.Boot = apply(LOADS, 2, mean),
	        Std.Error = apply(LOADS, 2, sd), perc.05 = apply(LOADS, 2, function(x) quantile(x, 0.05)),
	        perc.95 = apply(LOADS, 2, function(x) quantile(x, 0.95)))
	    colnames(PATHS) <- path.labs
	    PB <- data.frame(Original = path.orig, Mean.Boot = apply(PATHS, 2, mean),
	        Std.Error = apply(PATHS, 2, sd), perc.05 = apply(PATHS, 2, function(x) quantile(x, 0.05)),
	        perc.95 = apply(PATHS, 2, function(x) quantile(x, 0.95)))
	    colnames(TOEFS) <- Path.efs[, 1]
	    TE <- data.frame(Original = Path.efs[, 4], Mean.Boot = apply(TOEFS, 2, mean), 
	        Std.Error = apply(TOEFS, 2, sd), perc.05 = apply(TOEFS, 2, function(x) quantile(x, 0.05)), 
	        perc.95 = apply(TOEFS, 2, function(x) quantile(x, 0.95)))
	    colnames(RSQRS) <- lvs.names[endo == 1]
	    RB <- data.frame(Original = r2.orig, Mean.Boot = apply(RSQRS, 2, mean),
	        Std.Error = apply(RSQRS, 2, sd), perc.05 = apply(RSQRS, 2, function(x) quantile(x, 0.05)),
	        perc.95 = apply(RSQRS, 2, function(x) quantile(x, 0.95)))
	
	    # Bootstrap Results
	    res.boot[[signChangeCorrections[i]]] <- list(weights=WB, loadings=LB, paths=PB, rsq=RB, total.efs=TE)
	}

    return(res.boot)
}

.pls.efects <-
function(Path)
{
    lvs <- nrow(Path)
    lvs.names <- rownames(Path)
    path.efects <- as.list(1:(lvs-1))
    path.efects[[1]] <- Path
    if (lvs == 2)
    {
        ind.paths <- matrix(c(0,0,0,0),2,2)
        total.paths <- Path
    }
    if (lvs > 2)
    {
        for (k in 2:(lvs-1))
            path.efects[[k]] <- path.efects[[k-1]] %*% Path
        ind.paths <- matrix(0, lvs, lvs)
        for (k in 2:length(path.efects))
            ind.paths <- ind.paths + path.efects[[k]]
        total.paths <- Path + ind.paths
    }
    efs.labs <- NULL
    dir.efs <- NULL
    ind.efs <- NULL
    tot.efs <- NULL
    for (j in 1:lvs)
        for (i in j:lvs)
             if (i != j) 
             {
                 efs.labs <- c(efs.labs, paste(lvs.names[j],"->",lvs.names[i],sep=""))
                 dir.efs <- c(dir.efs, Path[i,j])# direct effects
                 ind.efs <- c(ind.efs, ind.paths[i,j])# indirect effects
                 tot.efs <- c(tot.efs, total.paths[i,j])# total effects
             }
    Effects <- data.frame(relationships=efs.labs, dir.effects=dir.efs, 
                          ind.effects=ind.efs, tot.effects=tot.efs)
    return(Effects)
}

.pls.gof <-
function(comu, R2, blocks, IDM)
{
    lvs <- nrow(IDM)
    blocklist <- as.list(1:lvs)
    for (j in 1:lvs)         
         blocklist[[j]] <- rep(j,blocks[j])
    blocklist <- unlist(blocklist)
    endo <- rowSums(IDM)
    endo[endo!=0] <- 1  
    n.end <- sum(endo)
    # average of communalities
    R2.aux <- R2[endo==1]
    comu.aux <- 0 
    n.comu <- 0
    for (j in 1:lvs)
    {
        if (length(which(blocklist==j))>1)
        {
            comu.aux <- comu.aux + sum(comu[which(blocklist==j)])
            n.comu <- n.comu + length(which(blocklist==j))
        }
    }
    gof <- sqrt((comu.aux/n.comu)*mean(R2.aux))
    return(gof)
}

.pls.GOF <-
function(DM, IDM, blocks, comu, unidim, R2)
{
    lvs <- nrow(IDM)
    blocklist <- as.list(1:lvs)
    for (j in 1:lvs)
         blocklist[[j]] <- rep(j,blocks[j])
    blocklist <- unlist(blocklist)
    endo <- rowSums(IDM)
    endo[endo!=0] <- 1  
    R2.aux <- R2[endo==1]
    comu.aux <- 0# storing communalities
    n.comu <- 0
    gof.out <- 0# storing outer model term
    gof.inn <- rep(1, lvs)  # inner model term
    for (j in 1:lvs)
    {
        k <- which(blocklist==j)
        if (length(k)>1)
        {
            comu.aux <- comu.aux + sum(comu[k])
            n.comu <- n.comu + length(k)
            gof.out <- gof.out + ( sum(comu[k]) / unidim[j,5] )
        }
        if (endo[j]==1)
        {
             EB <- DM[,which(blocklist==j)]  # endog block
             aux <- which(IDM[j,]==1)
             SB <- DM[,which(blocklist %in% aux)]  # super block             
             gof.inn[j] <- cancor(EB, SB)$cor[1]^2
        }
    }
    # ============================ GoF Indexes ===========================
    gof.abs <- sqrt((comu.aux/n.comu)*mean(R2.aux))
    gof.om <- sqrt(gof.out/length(which(blocks>1)))
    gof.im <- sqrt(sum(R2/gof.inn)/sum(endo))
    gof.rel <- gof.om * gof.im
    GOF <- data.frame(GoF=c("Absolute", "Relative", "Outer.mod", "Inner.mod"),
                      value=c(gof.abs, gof.rel, gof.om, gof.im))
    return(GOF)
}

.pls.GQI <-
function(pls, part, DM)
{
    # ========================== GQI function ==========================
    # Function to calculate Group Quality Index (GQI)        
    # =========================== arguments ==============================
    # pls: object of class "plspm"
    # part: vector with units memberships / or categorical variable
    # DM: data matrix
    
    IDM <- pls$model$IDM# Inner Design Matrix
    blocks <- pls$model$blocks# cardinality of blocks
    scheme <- pls$model$scheme# inner weighting scheme
    modes <- pls$model$modes# measurement modes
    scaled <- pls$model$scaled# type of scaling
    plsr <- FALSE 
    tol <- pls$model$tol
    iter <- pls$model$iter
    outer <- pls$model$outer
    blocklist <- outer
    for (k in 1:length(blocks))
         blocklist[[k]] <- rep(k,blocks[k])
    blocklist <- unlist(blocklist)
    lvs <- nrow(IDM)
    lvs.names <- rownames(IDM)
    mvs <- sum(blocks)
    blocklist <- as.list(1:lvs)
    for (j in 1:lvs)
         blocklist[[j]] <- rep(j,blocks[j])
    blocklist <- unlist(blocklist)
    endo <- rowSums(IDM)
    endo[endo!=0] <- 1  
    # data scaling (standardized data)
    sd.X <- sqrt((nrow(DM)-1)/nrow(DM)) * apply(DM, 2, sd)
    X <- scale(DM, scale=sd.X)
    clas.part <- part 
    # number of clusters
    nclus <- nlevels(factor(clas.part))
    w.locals <- as.list(1:nclus)# outer.weights
    LV.locals <- as.list(1:nclus)# std latent variables 
    loads.locals <- as.list(1:nclus)# loadings
    path.locals <- as.list(1:nclus)# path coefficients
    R2.locals <- as.list(1:nclus)# R2
    comu.locals <- as.list(1:nclus)# mvs communalities
    outres.locals <- as.list(1:nclus)# communality residuals
    innres.locals <- as.list(1:nclus)# structural residuals
    out.term <- as.list(1:nclus)# outer term for GQI
    inn.term <- as.list(1:nclus)# inner term for GQI
    gqi.locals <- rep(0,nclus)# pseudo-gqi for each class

    # define MV matrix for each initial class
    split.DM <- as.list(1:nclus)
    split.X <- as.list(1:nclus)
    for (k in 1:nclus)
        split.DM[[k]] <- DM[clas.part==k,]            

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
        LV.locals[[k]] <- split.X[[k]] %*% out.ws[[2]]
        # calculating path coefficients for each class
        pathmod <- .pls.paths(IDM, LV.locals[[k]], plsr)
        path.locals[[k]] <- pathmod[[2]]
        R2.locals[[k]] <- pathmod[[3]][endo==1]
        # calculating loadings and communalities for each class
        loadcomu <- .pls.loads(split.X[[k]], LV.locals[[k]], blocks)    
        loads.locals[[k]] <- loadcomu[[1]]
        comu.locals[[k]] <- loadcomu[[2]]
        # computation of communality residuals (squared)
        out.res <- split.X[[k]]
        for (j in 1:lvs)
        {
            q <- which(blocklist==j) 
            X.hat <- LV.locals[[k]][,j] %*% t(loads.locals[[k]][q])
            out.res[,q] <- (split.X[[k]][,q] - X.hat)^2# outer residuals
        }
        outres.locals[[k]] <- out.res
        # computation of inner residuals (squared)
        if (sum(endo)!=1)
            Y.hat <- LV.locals[[k]] %*% t(path.locals[[k]][endo==1,])   
        if (sum(endo)==1)
            Y.hat <- LV.locals[[k]] %*% path.locals[[k]][endo==1,]        
        innres.locals[[k]] <- (LV.locals[[k]][,endo==1] - Y.hat)^2
        # outer and inner terms of GQI formula
        out.term[[k]] <- mean((1-(colSums(outres.locals[[k]])/colSums(split.X[[k]]^2))))
        if (sum(endo)==1)
            inn.term[[k]] <- mean((1-(colSums(innres.locals[[k]])/sum(LV.locals[[k]][,endo==1]^2))))
        if (sum(endo)!=1)
            inn.term[[k]] <- mean((1-(colSums(innres.locals[[k]])/colSums(LV.locals[[k]][,endo==1]^2))))
        gqi.locals[k] <- out.term[[k]] * inn.term[[k]]
    }
    unit.prop <- unlist(lapply(split.X, nrow))/nrow(DM)# proportion of units in each class
    GQI <- sqrt(sum(gqi.locals * unit.prop))
    return(GQI)
}

.pls.loads <-
function(X, Y.lvs, blocks)
{
    lvs <- length(blocks)
    mvs <- ncol(X)
    blocklist <- as.list(1:lvs)
    for (j in 1:lvs)
         blocklist[[j]] <- rep(j,blocks[j])
    blocklist <- unlist(blocklist)
    loads <- rep(NA, mvs)
    comu <- rep(NA, mvs)
    for (j in 1:lvs)
        loads[blocklist==j] <- cor(X[,blocklist==j], Y.lvs[,j])
    comu <- loads^2
    names(loads) <- colnames(X)  
    names(comu) <- colnames(X)
    res.loads <- list(loads, comu)
    return(res.loads)
}

.pls.locals.test <-
function(X, pls, g)
{
    # =========================== ARGUMENTS ==============================
    # X: data matrix related with g
    # pls: an object of class "plspm"
    # g: a factor with 2 levels indicating the groups to be compared

    # ========================== INPUTS SETTING ==========================
    IDM <- pls$model$IDM# Inner Design Matrix
    blocks <- pls$model$blocks# cardinality of blocks
    scheme <- pls$model$scheme# inner weighting scheme
    modes <- pls$model$modes# measurement modes
    scaled <- pls$model$scaled# type of scaling
    plsr <- pls$model$plsr# pls-regression
    tol <- pls$model$tol# tolerance criterion
    iter <- pls$model$iter# max num iterations
    lvs <- nrow(IDM)
    lvs.names <- rownames(IDM)
    mvs <- sum(blocks)
    blocklist <- as.list(1:lvs)
    for (j in 1:lvs)
         blocklist[[j]] <- rep(j,blocks[j])
    blocklist <- unlist(blocklist)
    reps <- 100
    path.labs <- NULL
    efs.labs <- NULL
    for (j in 1:lvs)
        for (i in j:lvs)
             if (IDM[i,j]==1) 
                 path.labs <- c(path.labs, paste(lvs.names[j],"->",lvs.names[i],sep=""))

    # ====================== Group1 model estimation =====================
    g1.lab <- levels(g)[1]
    group1 <- which(g==levels(g)[1])
    ng1 <- length(group1)
    if(scaled) {
        sd.g1 <- sqrt((ng1-1)/ng1) * apply(X[group1,], 2, sd)
        X.g1 <- scale(X[group1,], scale=sd.g1) 
    } else {
        X.g1 <- scale(X[group1,], scale=FALSE)
    }
    wgs.g1 <- .pls.weights(X.g1, IDM, blocks, modes, scheme, tol, iter)
    if (is.null(wgs.g1)) stop("The algorithm is non convergent") 
    cor.XY <- cor(X.g1, X.g1%*%wgs.g1[[2]])
    w.sig <- rep(NA,lvs)
    for (k in 1:lvs) 
         w.sig[k] <- ifelse(sum(sign(cor.XY[which(blocklist==k),k]))<=0,-1,1)
    Y1.lvs <- X.g1 %*% wgs.g1[[2]] %*% diag(w.sig,lvs,lvs)
    dimnames(Y1.lvs) <- list(rownames(X.g1), lvs.names)
    # Path coefficients 
    pathmod.g1 <- .pls.paths(IDM, Y1.lvs, plsr)
    innmod.g1 <- pathmod.g1[[1]]
    Path.g1 <- pathmod.g1[[2]]
    R2.g1 <- pathmod.g1[[3]]    
    path.g1 <- as.vector(Path.g1[which(IDM==1)])
    names(path.g1) <- path.labs
    # calculating loadings and communalities for each class
    loadcomu <- .pls.loads(X.g1, Y1.lvs, blocks)   
    load.g1 <- loadcomu[[1]]
    # gof
    gof.g1 <- .pls.gof(load.g1^2, R2.g1, blocks, IDM)

    # ====================== Group2 model estimation =====================
    g2.lab <- levels(g)[2]
    group2 <- which(g==levels(g)[2])
    ng2 <- length(group2)
    if(scaled) {
        sd.g2 <- sqrt((ng2-1)/ng2) * apply(X[group2,], 2, sd)
        X.g2 <- scale(X[group2,], scale=sd.g2) 
    } else {
        X.g2 <- scale(X[group2,], scale=FALSE)
    }
    wgs.g2 <- .pls.weights(X.g2, IDM, blocks, modes, scheme, tol, iter)
    if (is.null(wgs.g2)) stop("The algorithm is non convergent") 
    cor.XY <- cor(X[group2,], X.g2%*%wgs.g2[[2]])
    w.sig <- rep(NA,lvs)
    for (k in 1:lvs) 
         w.sig[k] <- ifelse(sum(sign(cor.XY[which(blocklist==k),k]))<=0,-1,1)
    Y2.lvs <- X.g2 %*% wgs.g2[[2]] %*% diag(w.sig,lvs,lvs)
    dimnames(Y2.lvs) <- list(rownames(X.g2), lvs.names)
    # Path coefficients 
    pathmod.g2 <- .pls.paths(IDM, Y2.lvs, plsr)
    innmod.g2 <- pathmod.g2[[1]]
    Path.g2 <- pathmod.g2[[2]]
    R2.g2 <- pathmod.g2[[3]]    
    path.g2 <- as.vector(Path.g2[which(IDM==1)])
    names(path.g2) <- path.labs
    # calculating loadings and communalities for each class
    loadcomu <- .pls.loads(X.g2, Y2.lvs, blocks)   
    load.g2 <- loadcomu[[1]]
    # gof
    gof.g2 <- .pls.gof(load.g2^2, R2.g2, blocks, IDM)

    # ====================== Group Comparison =====================
    difpath.orig <- abs(path.g1-path.g2)
    difload.orig <- abs(load.g1-load.g2)
    difgof.orig <- abs(gof.g1 - gof.g2)
    group1 <- which(g==levels(g)[1])
    group2 <- which(g==levels(g)[2])
    ng1 <- length(group1)
    ng2 <- length(group2)
    difpath.perm <- matrix(0, reps, sum(IDM))
    difload.perm <- matrix(0, reps, mvs)
    difgof.perm <- rep(0, reps)
    for (i in 1:reps)# multigroup permutation
    {
        permu <- sample(1:(ng1+ng2), ng1+ng2)
        samg1 <- permu[1:ng1]
        samg2 <- permu[(ng1+1):(ng1+ng2)]
        if(scaled) {
            sd.g1 <- sqrt((ng1-1)/ng1) * apply(X[samg1,], 2, sd)
            sd.g2 <- sqrt((ng2-1)/ng2) * apply(X[samg2,], 2, sd)
            X.g1 <- scale(X[samg1,], scale=sd.g1) 
            X.g2 <- scale(X[samg2,], scale=sd.g2) 
        } else {
            X.g1 <- scale(X[samg1,], scale=FALSE)
            X.g2 <- scale(X[samg2,], scale=FALSE)
        }
        wgs.g1 <- .pls.weights(X.g1, IDM, blocks, modes, scheme, tol, iter)
        wgs.g2 <- .pls.weights(X.g2, IDM, blocks, modes, scheme, tol, iter)
        if (is.null(wgs.g1)) stop("Non convergence in bootstrap samples") 
        if (is.null(wgs.g2)) stop("Non convergence in bootstrap samples") 
        cor.XY <- cor(X.g1, X.g1%*%wgs.g1[[2]])
        w.sig <- rep(NA,lvs)
        for (j in 1:lvs) 
             w.sig[j] <- ifelse(sum(sign(cor.XY[blocklist==j,j]))<=0,-1,1)
        Y1.lvs <- X.g1 %*% wgs.g1[[2]] %*% diag(w.sig,lvs,lvs)
        cor.XY <- cor(X.g2, X.g2%*%wgs.g2[[2]])
        w.sig <- rep(NA,lvs)
        for (j in 1:lvs) 
             w.sig[j] <- ifelse(sum(sign(cor.XY[blocklist==j,j]))<=0,-1,1)
        Y2.lvs <- X.g2 %*% wgs.g2[[2]] %*% diag(w.sig,lvs,lvs)
        pathmod.g1 <- .pls.paths(IDM, Y1.lvs, plsr)
        paths.g1 <- pathmod.g1[[2]]    
        pathmod.g2 <- .pls.paths(IDM, Y2.lvs, plsr)
        paths.g2 <- pathmod.g2[[2]]
        loadcomu <- .pls.loads(X.g1, Y1.lvs, blocks)   
        loads.g1 <- loadcomu[[1]]
        loadcomu <- .pls.loads(X.g2, Y2.lvs, blocks)   
        loads.g2 <- loadcomu[[1]]  
        gofs.g1 <- .pls.gof(loads.g1^2, R2.g1, blocks, IDM)  
        gofs.g2 <- .pls.gof(loads.g2^2, R2.g2, blocks, IDM)
        # difference between groups
        pp1 <- as.vector(paths.g1[which(IDM==1)])
        pp2 <- as.vector(paths.g2[which(IDM==1)])
        difpath.perm[i,] <- abs(pp1 - pp2)
        difload.perm[i,] <- abs(loads.g1 - loads.g2)
        difgof.perm[i] <- abs(gofs.g1 - gofs.g2)
    }   
    # p-value for path coefficients
    path.perm <- difpath.orig 
    for (j in 1:sum(IDM))         
        path.perm[j] <- length(which(difpath.orig[j]<difpath.perm[,j])) + 1
    path.val <- (1/(reps+1))*path.perm 
    signi.path <- rep("no",length(path.val))
    signi.path[path.val<0.05] <- "yes"
    res.path <- round(cbind(path.g1, path.g2, difpath.orig, path.val), 4)
    res1 <- data.frame(res.path, signi.path)
    colnames(res1) <- c(paste(rep("Class",2),levels(g),sep="."), 
                         "diff.abs", "p.value", "sig.05")  
    # p-values for loadings
    load.perm <- difload.orig 
    for (j in 1:mvs)
        load.perm[j] <- length(which(difload.orig[j]<difload.perm[,j])) + 1
    load.val <- (1/(reps+1))*load.perm 
    signi.load <- rep("no",length(load.val))
    signi.load[load.val<0.05] <- "yes"
    res.load <- round(cbind(load.g1, load.g2, difload.orig, load.val), 4)
    res2 <- data.frame(res.load, signi.load)
    colnames(res2) <- c(paste(rep("Class",2),levels(g),sep="."), 
                         "diff.abs", "p.value", "sig.05")  
    # p-values for gof
    gof.perm <- length(which(difgof.orig<difgof.perm)) + 1
    gof.val <- (1/(reps+1))*gof.perm 
    signi.gof <- rep("no",length(gof.val))
    signi.gof[gof.val<0.05] <- "yes"
    res3 <- data.frame(round(gof.g1,4), round(gof.g2,4), 
               round(difgof.orig,4), round(gof.val,4), signi.gof)
    names(res3) <- c(paste(rep("Class",2),levels(g),sep="."), 
                         "diff.abs", "p.value", "sig.05")  
    # list with results 
    resul <- list(paths=res1, loadings=res2, gof=res3)
    return(resul)
}

.pls.paths <-
function(IDM, Y.lvs, plsr)
{
    lvs.names <- colnames(IDM)
    endo = rowSums(IDM)
    endo[endo!=0] <- 1  # vector indicating endogenous LVs
    innmod <- as.list(1:sum(endo))
    Path <- IDM
    residuals <- as.list(1:sum(endo))
    R2 <- rep(0,nrow(IDM))
    for (aux in 1:sum(endo)) 
    {
        k1 <- which(endo==1)[aux]    # index for endo LV
        k2 <- which(IDM[k1,]==1)     # index for indep LVs
        if (length(k2)>1 & plsr) {               
            path.lm <- .plsr1(Y.lvs[,k2], Y.lvs[,k1], nc=2)
            Path[k1,k2] <- path.lm$coeffs
            residuals[[aux]] <- path.lm$resid
            R2[k1] <- path.lm$R2[1]
            inn.val <- c(path.lm$R2[1], path.lm$cte, path.lm$coeffs)
            inn.lab <- c("R2", "Intercept", paste(rep("path_",length(k2)),names(k2),sep=""))
            names(inn.val) <- NULL
            innmod[[aux]] <- data.frame(concept=inn.lab, value=round(inn.val,4))
        }
        if (length(k2)==1 | !plsr) {
            path.lm <- summary(lm(Y.lvs[,k1] ~ Y.lvs[,k2]))
            Path[k1,k2] <- path.lm$coef[-1,1]
            residuals[[aux]] <- path.lm$residuals  
            R2[k1] <- path.lm$r.squared
            inn.val <- c(path.lm$r.squared, path.lm$coef[,1])
            inn.lab <- c("R2", "Intercept", paste(rep("path_",length(k2)),names(k2),sep=""))
            names(inn.val) <- NULL
            innmod[[aux]] <- data.frame(concept=inn.lab, value=round(inn.val,4))
        }
    }
    names(innmod) <- lvs.names[endo!=0]  
    names(R2) <- lvs.names
    res.paths <- list(innmod, Path, R2, residuals)
    return(res.paths)
}

.pls.unidim <-
function(DM, blocks, modes)
{
    lvs <- length(blocks) 
    lvs.names <- names(blocks)
    blocklist <- as.list(1:lvs)
    for (j in 1:lvs)
         blocklist[[j]] <- rep(j,blocks[j])
    blocklist <- unlist(blocklist)
    Mode <- modes
    Mode[modes=="A"] <- "Reflective"
    Mode[modes=="B"] <- "Formative"   
    obs <- nrow(DM)
    sdvf <- sqrt((nrow(DM)-1)/nrow(DM)) 
    # Unidimensionality
    Alpha <- rep(1, lvs)# Cronbach's Alpha for each block
    Rho <- rep(1, lvs)# D.G. Rho for each block
    eig.1st <- rep(1,lvs)# first eigenvalue
    eig.2nd <- rep(0,lvs)# second eigenvalue
    for (aux in 1:lvs) 
    {      
        if (blocks[aux] != 1) 
        { 
            # scaling data
            DM.block <- DM[,which(blocklist==aux)]
            stdev.X <- apply(DM.block, 2, sd) * sdvf 
            X.uni <- scale(DM.block, scale=stdev.X)
            if (nrow(X.uni)<ncol(X.uni)) {   # more columns than rows
                acp <- princomp(t(X.uni)) 
                X.rho <- t(X.uni)
            } else {   # more rows than columns
                acp <- princomp(X.uni)
                X.rho <- X.uni
            }
            if (modes[aux]=="A") 
            {
                p = ncol(X.uni)
                a.denom <- var(rowSums(X.uni)) * sdvf^2
                a.numer <- 2*sum(cor(X.uni)[lower.tri(cor(X.uni))])
                alpha <- (a.numer / a.denom) * (p/(p-1))
                Alpha[aux] <- ifelse(alpha<0,0,alpha)
                numer.rho <- colSums(cor(X.rho, acp$scores[,1]))^2
                denom.rho <- numer.rho + (p - colSums(cor(X.rho, acp$scores[,1])^2) )
                Rho[aux] <- numer.rho / denom.rho
            } else {  # modes[aux]=="B"
                Alpha[aux] <- 0
                Rho[aux] <- 0
            }
            eig.1st[aux] <- acp$sdev[1]^2
            eig.2nd[aux] <- acp$sdev[2]^2
        }
    }
    unidim <- data.frame(Type.measure=Mode, MVs=blocks, C.alpha=Alpha, 
                         DG.rho=Rho, eig.1st, eig.2nd)
    rownames(unidim) <- lvs.names
    return(unidim)
}

.pls.weights <-
function(X, IDM, blocks, modes, scheme, tol, iter, W=NULL)
{
    lvs <- nrow(IDM)
    mvs <- ncol(X)
    sdv <- sqrt((nrow(X)-1)/nrow(X))   # std.dev factor correction
    blocklist <- as.list(1:lvs)
    for (j in 1:lvs)
         blocklist[[j]] <- rep(j,blocks[j])
    blocklist <- unlist(blocklist)
    # outer design matrix 'ODM' and matrix of outer weights 'W'
    ODM <- matrix(0, mvs, lvs)
    for (j in 1:lvs)
        ODM[which(blocklist==j),j] <- rep(1,blocks[j])
     
    # If starting values of W are not specified, set them.
    
    if(is.null(W))
	    W <- ODM %*% diag(1/(sd(X %*% ODM)*sdv),lvs,lvs)
    w.old <- rowSums(W)    
    w.dif <- 1
    itermax <- 1
    repeat 
    {            
        Y <- X %*% W  # external estimation of LVs 'Y'
        Y <- scale(Y) * sdv
        # matrix of inner weights 'e' 
        E <- switch(scheme, 
               "centroid" = sign(cor(Y) * (IDM + t(IDM))),
               "factor" = cor(Y) * (IDM + t(IDM)),
               "path" = .path.scheme(IDM, Y))
        Z <- Y %*% E  # internal estimation of LVs 'Z'
        # scaling Z
        Z <- Z %*% diag(1/(sd(Z)*sdv), lvs, lvs)
        # computing outer weights 'w'
        for (j in 1:lvs)
        {
            X.blok = X[,which(blocklist==j)] 
            if (modes[j]=="A")# reflective way
                ODM[which(blocklist==j),j] <- (1/nrow(X)) * Z[,j] %*% X.blok
            if (modes[j]=="B")# formative way
                ODM[which(blocklist==j),j] <- solve.qr(qr(X.blok),Z[,j])
        }
        W <- ODM
        w.new <- rowSums(W)                
        w.dif <- sum((abs(w.old) - abs(w.new))^2)  # difference of out.weights 
        if (w.dif<tol || itermax==iter) break
        w.old <- w.new
        itermax <- itermax + 1
    } # end repeat       
    W <- ODM %*% diag(1/(sd(X %*% ODM)*sdv),lvs,lvs)
    w.new <- rowSums(W)                
    names(w.new) <- colnames(X)
    dimnames(W) <- list(colnames(X),rownames(IDM))       
    res.ws <- list(w.new, W, itermax)
    if (itermax==iter) res.ws=NULL
    return(res.ws)
}

.plsr1 <-
function(x, y, nc=NULL, scaled=TRUE)
{
    # ============ checking arguments ============
    X <- as.matrix(x)
    Y <- as.matrix(y)
    n <- nrow(X)
    p <- ncol(X)
    if (is.null(nc))
        nc <- p
    # ============ setting inputs ==============
    if (scaled) Xx<-scale(X) else Xx<-scale(X,scale=F)
    if (scaled) Yy<-scale(Y) else Yy<-scale(Y,scale=F)
    X.old <- Xx
    Y.old <- Yy
    Th <- matrix(NA, n, nc)# matrix of X-scores
    Ph <- matrix(NA, p, nc)# matrix of X-loadings
    Wh <- matrix(NA, p, nc)# matrix of raw-weights
    Uh <- matrix(NA, n, nc)# matrix of Y-scores
    ch <- rep(NA, nc)# vector of y-loadings
    # ============ pls regression algorithm ==============
    for (h in 1:nc)
    {
        w.old <- t(X.old) %*% Y.old / sum(Y.old^2)
        w.new <- w.old / sqrt(sum(w.old^2)) # normalization
        t.new <- X.old %*% w.new
        p.new <- t(X.old) %*% t.new / sum(t.new^2) 
        c.new <- t(Y.old) %*% t.new / sum(t.new^2)
        u.new <- Y.old / as.vector(c.new)
        Y.old <- Y.old - t.new%*%c.new# deflate y.old
        X.old <- X.old - (t.new %*% t(p.new))# deflate X.old
        Th[,h] <- t.new
        Ph[,h] <- p.new
        Wh[,h] <- w.new
        Uh[,h] <- u.new
        ch[h] <- c.new
    }
    Ws <- Wh %*% solve(t(Ph)%*%Wh)# modified weights
    Bs <- as.vector(Ws %*% ch) # std beta coeffs    
    Br <- Bs * (rep(sd(Y),p)/apply(X,2,sd))   # beta coeffs
    cte <- as.vector(mean(y) - Br%*%apply(X,2,mean))# intercept
    y.hat <- X%*%Br+cte# y predicted
    resid <- as.vector(Y - y.hat)# residuals
    R2 <- as.vector(cor(Th, Yy))^2  # R2 coefficients    
    names(Br) <- colnames(X)
    names(resid) <- rownames(Y)
    names(y.hat) <- rownames(Y)
    names(R2) <- paste(rep("t",nc),1:nc,sep="")
    res <- list(coeffs=Br, coef.std=Bs, cte=cte, R2=R2[1:nc], resid=resid, y.pred=y.hat)    
    return(res)
}

.rec.hclust <-
function(index, lwd=1, lty=1, col="black")
{
    # index: index of the current tree to draw
    members <- get('members', envir= ._a2r_envir) 
    bottom  <- get('bottom',  envir= ._a2r_envir) 
    if (index<0){ # it is a leaf
        if(is.null(members)){
           ._a2r_counter <<- ._a2r_counter + 1
           return(list(x=._a2r_counter, n=1))
        }
        else{
            cc <- ._a2r_counter
            mm <- members[-index]
            polygon(x=c(cc, cc+mm/2, cc+mm), y=c(bottom, 0, bottom),
                    col=col, border = col, lwd=lwd)
            ._a2r_counter <<- ._a2r_counter + mm
            return(list(x=cc+mm/2, n=mm))
        }
    }    
    h.m   <- ._a2r_hclu$height[index]
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~ do left
    index.l  <- ._a2r_hclu$merge[index,1]    
    h.l <- if(index.l<0) 0 else ._a2r_hclu$height[index.l]
    if (h.l<._a2r_height_cut & h.m > ._a2r_height_cut){
        ._a2r_group <<- ._a2r_group + 1
        col.l <- get("col.down",envir=._a2r_envir)[._a2r_group]
        lwd.l <- get("lwd.down",envir=._a2r_envir)
        lty.l <- get("lty.down",envir=._a2r_envir)
    }
    else{
        col.l <- col
        lwd.l <- lwd
        lty.l <- lty
    }
    out.l   <- .rec.hclust(index.l, col=col.l, lty=lty.l, lwd=lwd.l)
    x.l     <- out.l$x
    n.l     <- out.l$n
        
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~ do right
    index.r  <- ._a2r_hclu$merge[index,2]
    h.r <- if(index.r<0) 0 else ._a2r_hclu$height[index.r]
    if (h.r<._a2r_height_cut & h.m > ._a2r_height_cut){
        ._a2r_group <<- ._a2r_group + 1
        col.r <- get("col.down",envir=._a2r_envir)[._a2r_group]
        lwd.r <- get("lwd.down",envir=._a2r_envir)
        lty.r <- get("lty.down",envir=._a2r_envir)
    }
    else{
        col.r <- col
        lwd.r <- lwd
        lty.r <- lty
    }
    out.r   <- .rec.hclust(index.r, col=col.r, lty=lty.r, lwd=lwd.r)
    x.r     <- out.r$x
    n.r     <- out.r$n
        
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~ draw what you have to draw    
    type <- get("type",envir=._a2r_envir)
    x.m  <- (x.r + x.l) / 2  
    n    <- n.r + n.l
    x.b  <- (n.r * x.r + n.l * x.l) / n      
    knot.pos <- get("knot.pos",envir=._a2r_envir)     
    x <- switch(knot.pos, mean=x.m, left=x.l, right= x.r,
            random = x.l + runif(1)*(x.r-x.l), bary=x.b)
            
    if (type=="rectangle"){
        segments(x0  = c(x.l, x.l, x.r),
                 x1  = c(x.l, x.r, x.r),
                 y0  = c(h.l, h.m, h.r),
                 y1  = c(h.m, h.m, h.m),
                 col = col,
                 lty = lty,
                 lwd = lwd)
    }
    if (type =="triangle"){
        segments(x0  = c(x.l, x.r),
                 x1  = c(x  , x),
                 y0  = c(h.l, h.r),
                 y1  = c(h.m, h.m),
                 col = col,
                 lty = lty,
                 lwd = lwd)
    }
                        
    list(x=x, n=n)
}

.weightsplot <-
function(IDM, blocks, out.weights, arr.pos=arr.pos,
                  box.prop=box.prop, box.cex=box.cex, cex.txt=cex.txt, newdev=newdev)
{
    lvs <- nrow(IDM)
    ini.vec <- cumsum(blocks) - blocks + 1
    end.vec <- cumsum(blocks) 
    ## plot of outer weights
    for (k in 1:lvs)
    {
        num.mvs <- blocks[[k]]
        names.mvs <- names(out.weights)[ini.vec[k]:end.vec[k]]
        names.mvs <- c(names.mvs,rownames(IDM)[k])
        box.types <- c(rep("rect",num.mvs),"ellipse")
        MW <- matrix(0,num.mvs+1,num.mvs+1)
        MW.col <- MW
        MW[num.mvs+1,] <- c(out.weights[ini.vec[k]:end.vec[k]],0)
        MW.col[MW<0] <- "red"# negative out.weights in red
        MW.col[MW>0] <- "blue"   # positive out.weights in blue
        if (newdev)
            dev.new()
        plotmat(round(MW,4), curve=0, name=names.mvs, lwd=1, box.type=box.types, arr.width=0,
           arr.pos=arr.pos, arr.lcol=MW.col, arr.col=MW.col, box.prop=box.prop, box.cex=box.cex,  
           cex.txt=cex.txt, main=c(paste(rownames(IDM)[k]),"weights"))
    }
}

