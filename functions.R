#! /usr/bin/env Rscript

# This file contains all function definitions needed for the MapReduce job.

#library(Matrix)
#library(stringr)
#library(ggm)

# Needed for standardization and other convenience functions
library(QuantPsyc)

#library(monomvn)

# The PLS packages
library(plspm)
library(semPLS)

# library(Hmisc)

# Needed for generating data from covariance matrices
library(MASS)

#
# Generates a random PLS design as a lower triangular matrix.
#

generateRandomModel<-function(numberOfConstructs,expectedNumberOfOutgoingPaths){

	#Probability of path being specified is the expected number of paths divided by the number of possible paths
	pathProbability<-expectedNumberOfOutgoingPaths*numberOfConstructs/(numberOfConstructs^2-numberOfConstructs)
	
	#Create a vector of 1s and 0s to store the structural model
	populationModel <- matrix(runif(numberOfConstructs)<pathProbability,numberOfConstructs,numberOfConstructs)
	
	#Set the upper triangle and diagonal to NA
	populationModel[upper.tri(populationModel,diag=TRUE)] <- NA
	
	#Ensure that the model is a valid PLS design and return it
	
	return(ensureThatModelIsValid(populationModel))

}

#
# Ensures that the matrix is a valid PLS design. If not, then add paths until it is.
#

ensureThatModelIsValid<-function(model){	

	pathcount <- colSums(model,na.rm=TRUE) + rowSums(model,na.rm=TRUE)
	
	numberOfConstructs<-ncol(model)
	
	for(i in 1:numberOfConstructs){
	
		#All constructs must have at least one emitted or received path for the model to be valid.
		
		if(pathcount[i]==0){
	
			row<-sample(c(1:numberOfConstructs)[c(1:numberOfConstructs)!=i],1)
			if(row<i){
				col<-row
				row<-i
			}
			else{
				col<-i
			}
			model[row,col]<-1
		}
	}
	return(model)
}	

#
# Takes in a binary matrix specifying which paths are included and returns a population model
# with these paths set to random values from a given interval
#
# Returns a list of 
#

setPopulationModelPathValues <- function(populationModelWhichPaths,populationPathValues){
	
	
	continue<-TRUE
	print(populationModelWhichPaths)
	while(continue){
		populationModel<-populationModelWhichPaths*runif(populationModelWhichPaths^2,min=populationPathValues[1],max=populationPathValues[2])

		# Finally, we need to ensure that the model is valid. We do this by
		# generating the population covariance matrix for the models by applying  # tracing rules, then setting the variances to 1 and testing if this is
		# a valid correlation matrix
		# If this check fails, the population model path values must be
		# regenerated. (This is very rare, but can happen.)
	
		# Start with a diagonal matrix
		sigma<-diag(ncol(populationModel))
		
		# Fill in the lower diagonal with covariances.
		for(row in 2:nrow(populationModel)){
			for(col in 1:(row-1)){
				print(sigma)
				sigma[row,col]=sum(populationModel[row,]*sigma[,col],na.rm=TRUE)
			}
		}
		
		#Copy the lower diagonal to upper diagonal
		sigma <- sigma+ t(sigma) - diag(ncol(populationModel))
		
		#Test if this is a valid correlation matrix by checking that all the eigen values are positive
		
		continue<-min(eigen(sigma,only.value=TRUE)$values)<0
	}
	
	return(list(paths=populationModel,covariances=sigma))
}

#
# Generates a model that is tested by altering paths. The model is a 
# valid PLS design.
#

#TODO: This does not work correctly
# [5,]   NA   NA   NA   NA
# [6,]   NA   NA   NA    0
# [7,]   NA   NA   NA   NA
# [8,]    1    1   NA    1

generateTestedModel<-function(populationModelWhichPaths,omittedPathsShare,extraPaths){
	
	#Remove paths based on expected percentage of omitted paths

	testedModel<-populationModelWhichPaths*(runif(length(populationModelWhichPaths))>=omittedPathsShare)
	
	#Add paths based on expected number of extra paths.
	
	possibleNewPaths<-(length(testedModel)^2-length(testedModel))/2- sum(populationModelWhichPaths)
	expectedNewPaths<-extraPaths * (length(testedModel)-1)
	probabilityForNewPath <- expectedNewPaths/possibleNewPaths
	
	testedModel<-testedModel + ((!populationModelWhichPaths) & (runif(length(testedModel)) < probabilityForNewPath ))*lower.tri(populationModelWhichPaths)

	return(ensureThatModelIsValid(testedModel))
}

#
# Generates a data object that has elements "constructs", "indicators", and
# "factorLoadings"
#

generateData <- function(populationCovariances,sampleSize,indicatorCount,factorLoading,factorLoadingInterval,maxErrorCorrelation,methodVariance){

	constructCount<-nrow(populationCovariances)
	
	# Create the true values of the constructs.
	
	constructs <- mvrnorm(n=sampleSize,rep(0,constructCount),populationCovariances)

	#Create a method factor
	methodFactor<-rnorm(sampleSize,sd=sqrt(methodVariance))

	# Calculate indicators
	
factorLoadingValues<-factorLoading+(runif(constructCount*indicatorCount)-.5)*2*factorLoadingInterval

	# Calculate the model implied part of indicators
	
	indicatorBase<-NULL
	for(construct in 1:constructCount){
		for(indicator in 1:indicatorCount){
			indicatorBase<-cbind(indicatorBase,constructs[,construct]*factorLoadingValues[construct*(indicatorCount-1)+indicator])
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
	indicators<-Make.Z(indicators)
	
	names(indicators)<-paste("i",1:ncol(indicators),sep="")

	#T as True score
	names(constructs)<-paste("T",1:ncol(constructs),sep="")

	#Return the data
	
	return(list(constructs=as.data.frame(constructs), factorLoadings=factorLoadingValues, indicators=as.data.frame(indicators)))
}

#
# Estimates one a path model with regression and summed scales
#

estimateWithRegression<-function(model,data){
	
	constructCount=ncol(model)
	indicatorCount=ncol(data)/constructCount
	sampleSize=nrow(data)
	paths<-NULL
	
	#Form standardized summed scales
	
	sumscales<-data.frame(row.names =c(1:sampleSize))
	
	for ( i in 1:constructCount ){
		sumscales<-cbind(sumscales,data[c(((i-1)*indicatorCount+1):(i*indicatorCount))])
		
	}
	
	names(sumscales)<-paste("C",c(1:constructCount),sep="")

	#Use the generated scores to evaluate regression paths that were included in the model

	modelRowSums<-rowSums(model)
	
	for ( i in 1:constructCount ){

		#Only proceed with the regression if the variable is endogenous
		
		if(modelRowSums[i]>0){
			
			# Form the regression equation
			
			independents <-colnames(model)[which(as.logical(model[i,]))]
			dependent<-rownames(model)[i]
			
			formulastr<-paste(dependent,paste(independents,collapse=" + "),sep=" ~ ")
			
			formulaobj<-as.formula(formulastr)

			
			#Store summary of the results as object. We will next extract the coefficients and their standard errors from this object

			tempresults<-lm(formulaobj,data=sumscales)

			stdcoefficients<-summary(tempresults)$coefficients

			for(k in 2:nrow(stdcoefficients)){
				
				#"From","To","Value","StandardError","ModelingTechnique"
				newrow<-c(row.names(stdcoefficients)[k],dependent,stdcoefficients[[k,1]],stdcoefficients[[k,2]],"regress")

				paths<-rbind(paths,newrow)
			}
		}
	}
	
	#Return the construct scores and path estimates
	
	return(list(constructs=sumscales,paths=paths))
}

estimateWithSemPLS<-function(model,data){

	constructCount=ncol(model)
	indicatorCount=ncol(data)/constructCount
	sampleSize=nrow(data)
	paths<-NULL

	# Parameters that determine what kind of PLS run we do
	
	indicatorModes<-"A"
	weightingScheme<-"A"
	signCorrection<-"ConstructLevelChanges"
		
	#The models are matrixes in the form "from to"
	
	modelpaths<-c()
		
	for ( i in 1:constructCount ){
		for ( j in 1:constructCount){
			if(model[i,j]==1) {
				#Path from variable on the colum to variable on row
				modelpaths <-c(modelpaths,paste("C",j,sep=""),paste("C",i,sep=""))
			}
	    }
	}
	
	loadings<-c()
	
	for(i in 1:ncol(indicators)){
		loadings<-c(loadings,paste("C",ceiling(i/indicatorCount),sep=""),paste("i",i,sep=""))
	}

	#Inner model is a "from to" list
	
	innermodel<-matrix(modelpaths, ncol=2, nrow=length(modelpaths)/2,byrow = TRUE)

	outermodel <-matrix(data = loadings, ncol=2,nrow=length(loadings)/2,byrow = TRUE)

	semPLS <- sempls(plsm(indicators, strucmod=innermodel, measuremod=outermodel),indicators,maxit=100,E=weightingschemes)
			
	plspaths<-semPLS$coefficients[grep("beta",row.names(semPLS$coefficients)),]
			
	#Sort the list so that the paths are in alphabetical order
	plspaths<-plspaths[sort.list(row.names(plspaths)),]
			
	#Parse "To" and "From" as new colums for plspaths 
	plspaths$To<-sub(".*-> ","",plspaths$Path,perl=TRUE)
	plspaths$From<-sub(" ->.*","",plspaths$Path,perl=TRUE)

	#Store construct scores
			
	constructs<-semPLS$factor_scores
			
	# Then bootstrap the model.

	#Bootstrappping does not always converge, so we need to do exception handling
	
	semPLSboot<-NULL
				
	#Tolerate 5 errors and then move on. Most typical error is 10 consecutive convergence failures.

	counter<-0

	while(is.null(semPLSboot)&counter<=5){
		counter<-counter+1
		tryCatch(
			semPLSboot<-bootsempls(semPLS, method=signcorrections[signindex],nboot=100)
			,error = function(e){}
		)
	}	

	#Store the bootstrap results and PLS results

	if(!is.null(semPLSboot)){
	
		# Store the boot strapped estimates 
		
		boottable<-summary(semPLSboot)$table
	
		#Filter out everything but betas and sort
		boottable<-boottable[grep("beta",row.names(boottable)),]
				
		boottable<-boottable[sort.list(row.names(boottable)),]

		#"From","To","Value","StandardError","ModelingTechnique"
		paths<-cbind(plspaths$From,plspaths$To,boottable$Estimate,boottable$Std.Error, "sempls")
	}
	
	return(list(constructs=constructs,paths=paths))
}

estimateWithPlspm<-function(model,data){
	
	constructCount=ncol(model)
	indicatorCount=ncol(data)/constructCount
	
	#Start by constructing with inner and outer model matrixes. 
		
	outer <- list()
	
	for ( i in 1:constructCount ){
		outer[[i]] <- c(((i-1)*indicatorCount+1):(i*indicatorCount))
	}
	
	plsResults<-plspm(data,model,outer, rep("A",constructCount), scheme= "path", boot.val=TRUE)

	return(constructs=plsResults$latents,paths=plsResults$boot)
}