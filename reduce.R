#! /usr/bin/env Rscript
# The reducer receives a population model and tested model. It runs PLS and summed scales on all specified experimental conditions for these models and outputs the results as CSV. 

source("parameters.R")

options(warn=-1)

con <- file("stdin", open = "r")

#One line means population model- all tested models - one data specification combination. See prepare.R for details on the protocol.

while (length(line <- readLines(con, n = 1, warn = FALSE)) > 0) {
	
	# Parse the specification for this set of replications

	fullSpecification<-unlist(strsplit(line, split="\t"))
	
	# First 11 elements are scalars describing the parameters for the
	# replication. 
	
	counter<-specification[1]
	populationModelNumber<-specification[2]
	numberOfConstructsIndex<-specification[3]
	expectedNumberOfOutgoingPathsIndex<-specification[3]
	populationPathValuesIndex<-specification[5]
	sampleSizeIndex<-specification[6]
	indicatorCountIndexs<-specification[7]
	factorLoadingIndex<-specification[8]
	factorLoadingIntervalIndex<-specification[9]
	maxErrorCorrelationIndex<-specification[10]
	methodVarianceIndex<-specification[11]

	# The rest are matrices. These are population covariance, population paths,
	# and all tested models. Calculate the number and length of matrices and
	# parse the population model
	
	testedModelCount<-length(numberOfConstructs)*length(expectedNumberOfOutgoingPaths)*length(populationPathValues)
	matrixCount<-2+testedModelCount
	
	matrixLength<-(length(fullSpecification)-11)/matrixCount
	
	populationModel<-list(covariances=matrix(fullSpecification[12:(12+matrixLength)],ncol=sqrt(matrixLength)),paths=matrix(fullSpecification[(13+matrixLength):(12+matrixLength*2)],ncol=sqrt(matrixLength)))
	
	# Generate data for this set of replications
	
	data<-generateData(populationModel$covariances,sampleSizes[sampleSizeIndex],indicatorCounts[indicatorCountIndex],factorLoadings[factorLoadingIndex],factorLoadingIntervals[factorLoadingIntervalIndex],maxErrorCorrelations[maxErrorCorrelationIndex],methodVariances[methodVarianceIndex])
	
	# Test the models saving results
	
	results<-list()
	
	for(testedModelIndex in 1:testedModelCount){
		testedModel<-matrix(fullSpecification[(13+matrixLength*(testedModelIndex+1)):(12+matrixLength*(testedModelIndex+2))],ncol=sqrt(matrixLength))
		
		# Write code for testing the model
		
		results$regression <- estimateWithRegression(testedModel,data$indicators)
		
		results$plspm <- estimateWithPlspm(testedModel,data$indicators)
		
		results$sempls <- estimateWithSemPLS(testedModel,data$indicators)

	}
	
	# Calculate some statistics that need to be calculated over all analyses
	# done for this particular data set
	
	
	#Print the results. Start with the full specification string
	
	cat(fullspecification,"\n",sep="")
}

close(con)
