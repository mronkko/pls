#! /usr/bin/env Rscript

# This file contains all function definitions needed for the MapReduce job.

#library(semPLS)
#library(Matrix)
#library(stringr)
#library(ggm)
#library(QuantPsyc)
#library(monomvn)
#library(MASS)

#
# Generates a random PLS design as a lower triangular matrix.
#

generateRandomModel<-function(numberOfConstructs,expectedNumberOfOutgoingPaths,populationPathValues){

	#Probability of path being specified is the expected number of paths divided by the number of possible paths
	pathProbability<-expectedNumberOfOutgoingPaths*numberOfConstructs/(numberOfConstructs^2-numberOfConstructs)
	
	#Create a vector of 1s and 0s to store the structural model
	populationModel <- matrix(runif(numberOfConstructs)<pathProbability,numberOfConstructs,numberOfConstructs)
	
	#Set the upper triangle and diagonal to zeros
	populationModel <- populationModel*lower.tri(populationModel)
	
	#Ensure that the model is a valid PLS design
	populationModel<-ensureThatModelIsValid(model)
	
	return(populationModel)

}

#
# Ensures that the matrix is a valid PLS design. If not, then add paths until it is.
#

ensureThatModelIsValid<-function(model)	

	pathcount <- colSums(model) + rowSums(model)
	
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
	
	while(true){
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
				sigma[row,col]=populationModel[row,]*sigma[,col]
			}
		}
		
		#Copy the lower diagonal to upper diagonal
		sigma <- sigma+ t(sigma) - diag(ncol(populationModel))
		
		#Test if this is a valid correlation matrix by checking that all the eigen values are positive
		
		continue<-min(eigen(sigma,only.value=TRUE)$values)<0
	}
	
	return(list(paths=populationModel,covariances=sigma))
}