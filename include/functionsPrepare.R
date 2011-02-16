#! /usr/bin/env Rscript

# This file contains all function definitions needed for preparing the input file for MapReduce.

# Needed for generating data from covariance matrices
library(MASS)

#
# Generates a random PLS design as a lower triangular matrix.
#

generateRandomModel<-function(numberOfConstructs,expectedNumberOfOutgoingPaths){

	#Probability of path being specified is the expected number of paths divided by the number of possible paths
	pathProbability<-expectedNumberOfOutgoingPaths*numberOfConstructs/(numberOfConstructs^2-numberOfConstructs)
	
	#Create a vector of 1s and 0s to store the structural model
	populationModel <- matrix(runif(numberOfConstructs^2)<pathProbability,numberOfConstructs,numberOfConstructs)
	
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
	
	if(continue<-min(eigen(sigma,only.value=TRUE)$values)<0){
		return(NULL)
	}
	else{
		return(list(paths=populationModel,covariances=sigma))
	}
}

#
# Generates a model that is tested by altering paths. The model is a 
# valid PLS design.
#


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
