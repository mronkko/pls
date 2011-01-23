#! /usr/bin/env Rscript

# This file contains all function definitions needed for the MapReduce job.

#library(Matrix)
#library(stringr)
#library(ggm)


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
# Ensures that the matrix is a valid PLS design. If not, then add paths until it 
# is. There are two criterion that the model must fulfill
# 1) Each construct must be linked to at least one other with a regression path
# 2) It must be possible to go from one construct to every other by following 
# the regression paths
#

ensureThatModelIsValid<-function(model){	

	# Condition 1
	
	pathcount <- colSums(model,na.rm=TRUE) + rowSums(model,na.rm=TRUE)
	
	numberOfConstructs<-ncol(model)
	
	for(i in 1:numberOfConstructs){
	
		# All constructs must have at least one emitted or received path for the 
		# model to be valid. (Condition 1)
		
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
	
	# Condition 2. (There is probably a more elegent way to do this.)
	
	# Create a test matrix that records indirect relationships. If all 
	# constructs are connected, the all indirect relationships should be present 
	# if we convert each directional regression path to bidirectional.
	# If the model does not pass this test, add one non-redundant path and
	# repeat the test
	
	testMatrix<-model
	continue<-TRUE
	
	found<-FALSE
	while(continue){
		for(i in 1:numberOfConstructs){

			# Mark all variables that predict or depend on the i:th construct as 
			# connected in the testMatrix
		
			temp<-testMatrix[i,]|testMatrix[,i]
			testMatrix<-testMatrix | outer(temp,temp)
		}
		
		if(! is.na(min(testMatrix))){
			continue<-FALSE
		}
		else{
			# Choose one of the non-existing paths and add it to the model
			indexForNewPath<-sample(which(is.na(testMatrix)&lower.tri(testMatrix)),1)
			model[indexForNewPath]<-1
			testMatrix[indexForNewPath]<-TRUE
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
	
	possibleNewPaths<-(length(testedModel)-ncol(testedModel))/2- sum(populationModelWhichPaths,na.rm =TRUE)
	expectedNewPaths<-extraPaths * (ncol(testedModel)-1)
	probabilityForNewPath <- expectedNewPaths/possibleNewPaths

	
	testedModel<-testedModel + ((!populationModelWhichPaths) & (runif(length(testedModel)) < probabilityForNewPath ))*lower.tri(populationModelWhichPaths)

	testedModel<-ensureThatModelIsValid(testedModel)

	return(testedModel)
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
	indicators<-as.data.frame(scale(indicators))
	colnames(indicators)<-paste("i",1:ncol(indicators),sep="")

	#T as True score
	constructs<-as.data.frame(constructs)
	colnames(constructs)<-paste("T",1:ncol(constructs),sep="")

	#Return the data
	
	return(list(constructs=constructs, factorLoadings=factorLoadingValues, indicators=indicators))
}

#
# Estimates a path model with regression and summed scales
#

estimateWithRegression<-function(model,data){
	
	constructCount=ncol(model)
	indicatorCount=ncol(data)/constructCount
	sampleSize=nrow(data)
	paths<-NULL
	
	#Form standardized summed scales
	
	sumscales<-data.frame(row.names =c(1:sampleSize))
	
	for ( i in 1:constructCount ){
		sumscales<-cbind(sumscales,rowSums(data[,c(((i-1)*indicatorCount+1):(i*indicatorCount))]))
		
	}
	
	sumscales<-as.data.frame(scale(sumscales))
	colnames(sumscales)<-paste("C",c(1:constructCount),sep="")

	
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

			tempresults<-lm(formulaobj,data=sumscales)

			stdcoefficients<-summary(tempresults)$coefficients

			for(k in 2:nrow(stdcoefficients)){
				
				#"From","To","Estimates value","Mean.Boot","Std.Error","perc.05","perc.95","ModelingTechnique"
				newrow<-c(row.names(stdcoefficients)[k],dependent,stdcoefficients[[k,1]],NA,stdcoefficients[[k,2]],NA,NA,"regress")

				paths<-rbind(paths,newrow)
			}
		}
	}
	
	#Return the construct scores and path estimates
	
	return(list(constructs=sumscales,paths=paths))
}

#
# Estimates a model with SemPLS and bootstraps the model. This could probably be 
# written in much more compact form.
#


estimateWithSemPLS<-function(model,data){

	constructCount=ncol(model)
	indicatorCount=ncol(data)/constructCount
	paths<-NULL

		
	# The models are matrices in the form "from to" 
	
	modelpaths<-c()
		
	for ( i in 2:constructCount ){
		for ( j in 1:(i-1)){
			if(model[i,j]==1) {
				#Path from variable on the colum to variable on row
				modelpaths <-c(modelpaths,paste("C",j,sep=""),paste("C",i,sep=""))
			}
	    }
	}
	
	# Loadings are matrices in the from "from to"
	loadings<-c()
	
	for(i in 1:ncol(data)){
		loadings<-c(loadings,paste("C",ceiling(i/indicatorCount),sep=""),paste("i",i,sep=""))
	}

	# construct the matrices in the right format
	
	innermodel<-matrix(modelpaths, ncol=2, nrow=length(modelpaths)/2,byrow = TRUE)

	outermodel <-matrix(loadings, ncol=2,nrow=length(loadings)/2,byrow = TRUE)

	semPLS <- sempls(plsm(data, strucmod=innermodel, measuremod=outermodel),data,maxit=300)
			
	plspaths<-semPLS$coefficients[grep("beta",row.names(semPLS$coefficients)),]
			
	#Sort the list so that the paths are in alphabetical order
	plspaths<-plspaths[sort.list(row.names(plspaths)),]
			
	#Parse "To" and "From" as new colums for plspaths 
	plspaths$To<-sub(".*-> ","",plspaths$Path,perl=TRUE)
	plspaths$From<-sub(" ->.*","",plspaths$Path,perl=TRUE)

	#Store construct scores
			
	constructs<-semPLS$factor_scores
			
	# Then bootstrap the model.

	
	semPLSboot<-NULL

	# Bootstrappping does not always converge, so we need to do exception 
	# handling. Tolerate 5 errors and then move on. Most typical error is 10 
	# consecutive convergence failures.

	counter<-0

	while(is.null(semPLSboot)&counter<=5){
		counter<-counter+1
		tryCatch(
			semPLSboot<-bootsempls(semPLS,nboot=100)
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
	
	# plspm does not like NA:s in the matrix, so we will replace these with 
	# zeros.
	model[is.na(model)]<-0
	
	plsResults<-plspm(data,model,outer, rep("A",constructCount), scheme= "path", boot.val=TRUE)

	# "From","To","Estimated value","Mean.Boot","Std.Error","perc.05","perc.95","ModelingTechnique"
	paths<-cbind(sub("->.*","",rownames(plsResults$boot$paths)),sub(".*->","",rownames(plsResults$boot$paths)),plsResults$boot$paths,"pls")

	return(list(constructs=plsResults$latents,paths=paths))
}

#
# A wrapper for print to allow easily commenting out all unneccessary print commmands
#
debugPrint<-function(x){
	print(x)
}