#! /usr/bin/env Rscript
# The reducer receives a population model and tested model. It runs PLS and summed scales on all specified experimental conditions for these models and outputs the results as CSV. 

source("parameters.R")

# Some packages that we include want to write things to stdout. Redirect this to
# /dev/null

source("functions.R")

options(warn=-1)

con <- file("stdin", open = "r")

#One line means population model- all tested models - one data specification combination. See prepare.R for details on the protocol.

while (length(line <- readLines(con, n = 1, warn = FALSE)) > 0) {
	
	# Parse the specification for this set of replications

	specification<-sapply(unlist(strsplit(line, split="\t")),as.numeric)
	
	# Setting colnames will help in debugging
	names(specification)<-c(1:length(specification))
	
	# First 11 elements are scalars describing the parameters for the
	# replication. 
	
	counter<-specification[1]
	populationModelNumber<-specification[2]
	numberOfConstructsIndex<-specification[3]
	expectedNumberOfOutgoingPathsIndex<-specification[3]
	populationPathValuesIndex<-specification[5]
	sampleSizeIndex<-specification[6]
	indicatorCountIndex<-specification[7]
	factorLoadingIndex<-specification[8]
	factorLoadingIntervalIndex<-specification[9]
	maxErrorCorrelationIndex<-specification[10]
	methodVarianceIndex<-specification[11]

	constructCount<-numberOfConstructs[numberOfConstructsIndex]
	
	# The rest are matrices. These are population covariance, population paths,
	# and all tested models. Calculate the number and length of matrices and
	# parse the population model
	
	testedModelCount<-length(expectedNumberOfOutgoingPaths)*length(populationPathValues)
	matrixCount<-2+testedModelCount
	
	matrixLength<-constructCount^2

	populationModel<-list(covariances=matrix(specification[12:(11+matrixLength)],ncol=constructCount),paths=matrix(specification[(12+matrixLength):(11+matrixLength*2)],ncol=constructCount))
	
	# Generate data for this set of replications

	data<-generateData(populationModel$covariances,sampleSizes[sampleSizeIndex],indicatorCounts[indicatorCountIndex],factorLoadings[factorLoadingIndex],factorLoadingIntervals[factorLoadingIntervalIndex],maxErrorCorrelations[maxErrorCorrelationIndex],methodVariances[methodVarianceIndex])
	
	# Test the models saving results
	
	results<-list()
	
	print(specification)
	
	for(testedModelIndex in 1:testedModelCount){
		results[testedModelIndex]<-list()
		testedModel<-matrix(specification[(13+matrixLength*(testedModelIndex+1)):(12+matrixLength*(testedModelIndex+2))],ncol=sqrt(matrixLength))
		
		print(paste(testedModelIndex,13+matrixLength*(testedModelIndex+1),12+matrixLength*(testedModelIndex+2)))
		print(testedModel)
		
		# Model estimation requires that the names are correctly included
		colnames(testedModel)<-paste("C",c(1:constructCount),sep="")
		rownames(testedModel)<-paste("C",c(1:constructCount),sep="")
		
		# Write code for testing the model
		
		# results[testedModelIndex]$regression <- estimateWithRegression(testedModel,data$indicators)
		
		# results[testedModelIndex]$plspm <- estimateWithPlspm(testedModel,data$indicators)
		
		# results[testedModelIndex]$sempls <- estimateWithSemPLS(testedModel,data$indicators)

	}
	
	# Calculate some statistics that need to be calculated over all analyses
	# done for this particular data set
	
	#
	# Calculate correlations for each construct over all estimated models.
	# We do this by calculating a mean score over all nine tests for each type
	# of test and then test calculating how well each replication correlates 
	# with the mean values
	#

	for(modelTypeIndex in 1:length(results)){
	
		constructSums<-NULL
		
		for(testedModelIndex in 1:testedModelCount){
			constructSums<-constructSums+results[[testedModelIndex]][[modelTypeIndex]]$constructs
		}
		constructMeans[[modelTypeIndex]]=construcSums/length(results)
	}
	
	
	#Print the results. Start with the full specification string (input line)
	
	cat(line,"\n",sep="")

	#Print true factor loadings
	
	#Loop over replications and model types. Print the resutls as lines.
	
	for(modelTypeIndex in 1:length(results)){
	
		for(testedModelIndex in 1:testedModelCount){
			
			#Print a matrix of path coefficients
			
			#Print a matrix of standard errors
			
			#Print a correlation matrix of construct scores and true scores
			
			#Print item-construct correlation matrix
			
			#Print a vector of correlations with mean of construct estimates
		}
	}

}

close(con)
