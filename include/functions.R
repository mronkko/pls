#! /usr/bin/env Rscript

# This file contains all function definitions needed for the MapReduce job.

# The PLS packages. Instead of loading the library, load modified files from cache
# library(plspm)

source("include/plspm-internal.R")
source("include/plspm.R")

#
# A wrapper for print to allow easily commenting out all unneccessary print commmands
#

debugPrint<-function(x){
	print(x)
}

#
# mat.regress from package psych. Calculates regression from covariance matrix. We save a small 
# amount of computing time by copying this function instead of installing the psych library.
#
#

mat.regress <- function (m, x, y, n.obs = NULL) 
{
    cl <- match.call()
    if (!is.matrix(m)) 
        m <- as.matrix(m)
    if (dim(m)[1] != dim(m)[2]) {
        n.obs = dim(m)[1]
        C <- cov(m, use = "pairwise")
        m <- cov2cor(C)
    }
    else {
        C <- m
        m <- cov2cor(m)
    }
    nm <- dim(m)[1]
    xy <- c(x, y)
    numx <- length(x)
    numy <- length(y)
    nxy <- numx + numy
    a.matrix <- m[x, x]
    ac.matrix <- C[x, x]
    b.matrix <- m[x, y]
    bc.matrix <- C[x, y]
    if (numx == 1) {
        beta <- matrix(b.matrix, nrow = 1)
    }
    else {
        beta <- solve(a.matrix, b.matrix)
        beta <- as.matrix(beta)
    }
    if (numy > 1) {
        if (is.null(rownames(beta))) {
            rownames(beta) <- colnames(m)[x]
        }
        if (is.null(colnames(beta))) {
            colnames(beta) <- colnames(m)[y]
        }
        R2 <- colSums(beta * b.matrix)
    }
    else {
        colnames(beta) <- colnames(m)[1]
        R2 <- sum(beta * b.matrix)
        names(beta) <- colnames(m)[x]
        names(R2) <- colnames(m)[y]
    }
    if (!is.null(n.obs)) {
        k <- length(x)
        uniq <- (1 - smc(a.matrix))
        se.beta <- list()
        for (i in 1:length(y)) {
            df <- n.obs - k - 1
            se.beta[[i]] <- (sqrt((1 - R2[i])/(df)) * sqrt(1/uniq))
        }
        se <- matrix(unlist(se.beta), ncol = length(y))
        colnames(se) <- colnames(beta)
        rownames(se) <- rownames(beta)
        tvalue <- beta/se
        se <- t(t(se) * sqrt(diag(C)[y]))/sqrt(diag(ac.matrix))
        prob <- 2 * (1 - pt(abs(tvalue), df))
        SE2 <- 4 * R2 * (1 - R2)^2 * (df^2)/((n.obs^2 - 1) * 
            (n.obs + 3))
        SE = sqrt(SE2)
        F <- R2 * df/(k * (1 - R2))
        pF <- 1 - pf(F, k, df)
        shrunkenR2 <- 1 - (1 - R2) * (n.obs - 1)/df
    }
    if (numx == 1) {
        beta <- beta * sqrt(diag(C)[y])
    }
    else {
        beta <- t(t(beta) * sqrt(diag(C)[y]))/sqrt(diag(ac.matrix))
    }
    if (is.null(n.obs)) {
        mat.regress <- list(beta = beta, R = sqrt(R2), R2 = R2, 
            Call = cl)
    }
    else {
        mat.regress <- list(beta = beta, se = se, t = tvalue, 
            Probability = prob, R = sqrt(R2), R2 = R2, shrunkenR2 = shrunkenR2, 
            seR2 = SE, F = F, probF = pF, df = c(k, df), Call = cl)
    }
    class(mat.regress) <- c("psych", "mat.regress")
    return(mat.regress)
}

#
# Draws a sample from multivariate normal distribution with the given covariance matrix
#
# This is from MASS package, version 3.7-9. Since we need only this one function, it does not
# make sense to install the entire package, which is quite large.
#

mvrnorm <- function (n = 1, mu, Sigma, tol = 1e-06, empirical = FALSE) 
{
    p <- length(mu)
    if (!all(dim(Sigma) == c(p, p))) 
        stop("incompatible arguments")
    eS <- eigen(Sigma, symmetric = TRUE, EISPACK = TRUE)
    ev <- eS$values
    if (!all(ev >= -tol * abs(ev[1L]))) 
        stop("'Sigma' is not positive definite")
    X <- matrix(rnorm(p * n), n)
    if (empirical) {
        X <- scale(X, TRUE, FALSE)
        X <- X %*% svd(X, nu = 0)$v
        X <- scale(X, FALSE, TRUE)
    }
    X <- drop(mu) + eS$vectors %*% diag(sqrt(pmax(ev, 0)), p) %*% 
        t(X)
    nm <- names(mu)
    if (is.null(nm) && !is.null(dn <- dimnames(Sigma))) 
        nm <- dn[[1L]]
    dimnames(X) <- list(nm, NULL)
    if (n == 1) 
        drop(X)
    else t(X)
}


#
# Generates a data object that has elements "constructs", "indicators", and
# "factorLoadings"
#

generateData <- function(constructs,indicatorCount,factorLoading,factorLoadingInterval,maxErrorCorrelation,methodVariance){

	constructCount<-ncol(constructs)
	sampleSize<-nrow(constructs)
	
	#Create a method factor
	methodFactor<-rnorm(sampleSize,sd=sqrt(methodVariance))

	# Calculate indicators
	factorLoadingValues<-factorLoading+(runif(constructCount*indicatorCount)-.5)*2*factorLoadingInterval

	# Calculate the model implied part of indicators
	
	indicatorBase<-NULL
	for(construct in 1:constructCount){
		for(indicator in 1:indicatorCount){
			indicatorBase<-cbind(indicatorBase,constructs[,construct]*factorLoadingValues[(construct-1)*indicatorCount+indicator])
		}
	}
	
	# Calculate how much variance is left for error variance
	
	errorTermVariances <- 1-factorLoadingValues^2
	
	# Add method variances. If method variance would cause the total variance to # go over one in the population, scale down the method variance for the 
	# particular indicator.
	
	methodVariances=min(methodVariance,errorTermVariances)
	
	methodVarianceComponent <- matrix(rep(rnorm(sampleSize) %o% methodVariances,constructCount*indicatorCount),ncol=constructCount*indicatorCount)
	
	# Update error term variances now that we have added method variance

	errorTermVariances<- errorTermVariances-methodVariances

	# Then random errors, that are not nescessarily random because they can
	# correlate in the population

	# Start with random correlation matrix
	R <- matrix((runif((constructCount*indicatorCount)^2)*2-1)*.1, ncol=constructCount*indicatorCount)
	RtR <- R %*% t(R)
	R<-cov2cor(RtR)	
	# Scale the correlations with the maximum error correlation. Keep the diagonals as one
	
	R<-R*maxErrorCorrelation+(1-maxErrorCorrelation)*diag(indicatorCount*constructCount)
	
	# Generate a transformation matrix that scales the correlation matrix to be # a covariance matrix with the SDs of disturbances on the diagonal
	
	transmat<-sqrt(errorTermVariances) %*% t(sqrt(errorTermVariances))
	
	#Scale the correlation matrix to be error covariance matrix
	errorCovarianceMatrix <- R * transmat

	# Generate
	# error terms
	
	errorTerms <- mvrnorm(n=sampleSize,rep(0,constructCount*indicatorCount),errorCovarianceMatrix)

	
	# Indicators are asum of the components
	indicators=data.frame(indicatorBase+methodVarianceComponent+errorTerms)

	# Standardize the indicators
	indicators<-as.data.frame(scale(indicators))
	colnames(indicators)<-paste("i",1:ncol(indicators),sep="")

	#T as True score
	constructs<-as.data.frame(constructs)
	colnames(constructs)<-paste("T",1:ncol(constructs),sep="")

	#Return the data
	
	return(list(constructs=constructs, factorLoadings=factorLoadingValues, indicators=indicators))
}

#
# Estimates a path model with regression. Method is the name of the method used to generate the
# construct scores. (sumscale, factor, component)
#

estimateWithRegression<-function(model,data,method){
	
	constructCount=ncol(model)
	indicatorCount=ncol(data)/constructCount
	sampleSize=nrow(data)
	paths<-NULL
	
	#Form the construct scores
	
	
	constructScores<-data.frame(row.names =c(1:sampleSize))
	
		for ( i in 1:constructCount ){
			
			indicatorIndices<-c(((i-1)*indicatorCount+1):(i*indicatorCount))
			
			if(method=="sumscale"){
				constructScores<-cbind(constructScores,rowSums(data[,indicatorIndices]))
			}
			else if (method == "component"){
				
				# R tends to give principal components that have negative loadings.
				# Reverse the component so that the loadings are positive.
				
				obj<-princomp(data[,indicatorIndices],scores = TRUE)
				constructScores<-cbind(constructScores,obj$scores[,1]*((mean(obj$loadings[,1])>0)*2-1))
			}
			else if (method == "factor"){
				constructScores<-cbind(constructScores,factanal(data[,indicatorIndices],1,scores="regression")$scores[,1])
			}
	}
	
	# Standardize the scores and add names
	constructScores<-as.data.frame(scale(constructScores))
	colnames(constructScores)<-paste("C",c(1:constructCount),sep="")

	
	#Use the generated scores to evaluate regression paths that were included in the model

	modelRowSums<-rowSums(model,na.rm = TRUE)
	
	for ( i in 1:constructCount ){

		#Only proceed with the regression if the variable is endogenous
		
		if(modelRowSums[i]>0){
			
			# Form the regression equation
			
			independents <-colnames(model)[which(as.logical(model[i,]))]
			dependent<-rownames(model)[i]
			
			formulastr<-paste(dependent,paste(independents,collapse=" + "),sep=" ~ ")
			
			formulaobj<-as.formula(formulastr)

			
			#Store summary of the results as object. We will next extract the coefficients and their standard errors from this object

			tempresults<-lm(formulaobj,data=constructScores)

			stdcoefficients<-summary(tempresults)$coefficients

			for(k in 2:nrow(stdcoefficients)){
				
				#"From","To","Estimates value","Mean.Boot","Std.Error","perc.05","perc.95","ModelingTechnique"
				newrow<-c(row.names(stdcoefficients)[k],dependent,stdcoefficients[[k,1]],NA,stdcoefficients[[k,2]],NA,NA)

				paths<-rbind(paths,newrow)
			}
		}
	}
	
	colnames(paths)<-c("From","To","Estimate","Mean.Boot","Std.Error","perc.05","perc.95")
	rownames(paths)<-c(1:nrow(paths))

	#Return the construct scores and path estimates
	
	return(list(constructs=constructScores,paths=paths))
}

estimateWithPlspm<-function(model,data,doBootstrap){

	# Uncomment this to run wihtout PLS
	# return(NA)
	
	constructCount=ncol(model)
	indicatorCount=ncol(data)/constructCount
	
	#Start by constructing with inner and outer model matrixes. 
		
	outer <- list()
	
	for ( i in 1:constructCount ){
		outer[[i]] <- c(((i-1)*indicatorCount+1):(i*indicatorCount))
	}
	
	# plspm does not like NA:s in the matrix, so we will replace these with 
	# zeros.
	model[is.na(model)]<-0

	plsResults<-NULL
	
	plsResults<-plspm(data,model,outer, rep("A",constructCount), scheme= "path", iter=500, boot.val=doBootstrap)
	
	# "From","To","Estimated value","Mean.Boot","Std.Error","perc.05","perc.95","ModelingTechnique"
	
	if(doBootstrap){
		pathsStandard<-cbind(sub("->.*","",rownames(plsResults$boot$paths)),sub(".*->","",rownames(plsResults$boot$paths)),plsResults$boot.all[["Standard"]]$paths)
		pathsIndicatorCorrection<-cbind(sub("->.*","",rownames(plsResults$boot$paths)),sub(".*->","",rownames(plsResults$boot$paths)),plsResults$boot.all[["IndividualSignChanges"]]$paths)
		pathsConstructCorrection<-cbind(sub("->.*","",rownames(plsResults$boot$paths)),sub(".*->","",rownames(plsResults$boot$paths)),plsResults$boot.all[["ConstructLevelChanges"]]$paths)
	

		colnames(pathsStandard)<-c("From","To","Estimate","Mean.Boot","Std.Error","perc.05","perc.95")
		rownames(pathsStandard)<-c(1:nrow(pathsStandard))
		colnames(pathsIndicatorCorrection)<-c("From","To","Estimate","Mean.Boot","Std.Error","perc.05","perc.95")
		rownames(pathsIndicatorCorrection)<-c(1:nrow(pathsIndicatorCorrection))
		colnames(pathsConstructCorrection)<-c("From","To","Estimate","Mean.Boot","Std.Error","perc.05","perc.95")
		rownames(pathsConstructCorrection)<-c(1:nrow(pathsConstructCorrection))

		return(list(Standard=list(constructs=plsResults$latents,paths=pathsStandard),
				IndividualSignChanges=list(constructs=plsResults$latents,paths=pathsIndicatorCorrection),
				ConstructLevelChanges=list(constructs=plsResults$latents,paths=pathsConstructCorrection)))
	}

	else{
		pathsStandard<-NULL
		for(c in 1:ncol(plsResults$path.coefs)){
			for(r in 1:nrow(plsResults$path.coefs)){
				if(plsResults$path.coefs[r,c]!=0){
					pathsStandard<-rbind(pathsStandard,c(colnames(plsResults$path.coefs)[c],rownames(plsResults$path.coefs)[r],plsResults$path.coefs[r,c]))
				}
			}
		}
		colnames(pathsStandard)<-c("From","To","Estimate")
		rownames(pathsStandard)<-c(1:nrow(pathsStandard))

		return(list(Standard=list(constructs=plsResults$latents,paths=pathsStandard),IndividualSignChanges=NA,ConstructLevelChanges=NA))
	
	}
}


#
# This function generates the fractional factorial design matrix that is used in the simulation
# It is called fom the prepare.R and reducer.R so that we do not need to pass all simulation
# parameters to reducers, but only need to pass the design numbers. (Row indices for this matrix)
#

createdesignMatrix <- function(){
	
	
	# See the this article for details:
	# Xu, H. 2005. “A catalogue of three-level regular fractional factorial
	# designs,” Metrika (62:2), pp. 259–281.
	#
	# The key problem with creating the design is that the simulation runs are 
	# nested: We have several sets of construct scores. We want to run at least
	# three different measurements for each score while keeping the tested model
	# constant to estimate test-retest stability of the estimated (Constraint 1).
	# Also we need to use each indicator data to estimate at least three different 
	# models so that we can analyze the stability of the measurement (Constraint 2).
	#
	# All replications must be independent, so these two constraints need to be
	# taken care of in the design matrix. 
	#
	# We should have one identity column that does not affect the construt scores to # be able to estimate three models for each set of constucts. This identity 
	# column should not affect the indicator columns either to keep both the
	# indicator data constant (no changes in constructs or measurements) over three # model estimations
	#
	# The first design that fills this criteria is 11-5.10 from the 729 replications # set from the article above.
	#
	# See http://www.stat.ucla.edu/~hqxu/pub/ffd/ffd729.pdf
	#
	
	
	#
	# 729 runs would allow 14 factors with resolution V (14-8.1)
	
	generator <- matrix(c(	1,0,0,1,
							0,1,1,2),ncol=4,byrow=TRUE)
	
	for(i in 3:6){
		generator<-cbind(rbind(generator,0),c(rep(0,i-1),1),rbind(generator,1),rbind(generator,2))
		
	}
	
	
	#
	# This command prints the design matrix. The output is included below the 
	# command.
	#
	# print(generator[,c(1,2,5,14,41,122,63,149,166,188,222)])
	#
	#      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11]
	# [1,]    1    0    0    0    0    0    0    1    0     0     0
	# [2,]    0    1    0    0    0    0    1    2    1     1     0
	# [3,]    0    0    1    0    0    0    1    2    0     2     1
	# [4,]    0    0    0    1    0    0    1    1    0     1     1
	# [5,]    0    0    0    0    1    0    1    0    1     1     2
	# [6,]    0    0    0    0    0    1    0    1    1     1     1
	#
	# The factors are assigned as follows
	# Tested model: 1, 8
	# Measurement: 2,7,9,10,11
	# Population model and true scores: 3,4,5,6
	#

	# Get the indices from the 11-5.10 design reported in the online appendix of Xu,
	# 2005
	
	useGeneratorIndices<-c(1,2,5,14,41,122,63,149,166,188,222)
	generatorIndicesForPopulation<-c(2,3,4)
	generatorIndicesForSample<-c(6)
	generatorIndicesForTestedModel<-c(1,8)
	generatorIndicesForMeasurement<-c(5,7,9,10,11)
	
	#
	# permutations of a b c d e f
	#
	permutations<-matrix(c(c(0:728)%%3,c(0:728)%/%3%%3,c(0:728)%/%9%%3,c(0:728)%/%27%%3,c(0:728)%/%81%%3,c(0:728)%/%243),ncol=6)
	
	
	# Design matrix is the product of permutations and generator matrix, mod 3
	
	fullDesignMatrix<- ( permutations %*% generator ) %% 3
	
	
	designMatrix<-fullDesignMatrix[,useGeneratorIndices[c(	generatorIndicesForPopulation,
															generatorIndicesForSample,
															generatorIndicesForTestedModel,
															generatorIndicesForMeasurement)]]+1
	
	# Sort by identity columns and return

	designMatrix<-designMatrix[order(designMatrix[,1],designMatrix[,2],designMatrix[,3],designMatrix[,4],designMatrix[,5],designMatrix[,7]),]
	designMatrix<-cbind(designMatrix,matrix(c(1:729),ncol=1))
	
	colnames(designMatrix)<-c("numberOfConstructs","expectedNumberOfOutgoingPaths","populationPathValues","omittedPathsShare","extraPaths","sampleSize","indicatorCount","factorLoading","factorLoadingInterval","maxErrorCorrelation","methodVariance","designNumber")

	return(designMatrix)
	
}

#
# A small utility function to check memory use by object
#

displayMemory<-function(obs){
	for(i in 1:length(obs)){
		print(obs[i])
		print(object.size(get(obs[i])),units="Mb")
	}
}