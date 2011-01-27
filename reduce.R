#! /usr/bin/env Rscript
# The reducer receives a population model and tested model. It runs PLS and summed scales on all specified experimental conditions for these models and outputs the results as CSV. 

source("parameters.R")

source("functions.R")

options(warn=-1)

con <- file("stdin", open = "r")

#One line means population model- all tested models - one data specification combination. See prepare.R for details on the protocol.

while (length(line <- readLines(con, n = 1, warn = FALSE)) > 0) {

	# Record the time that this iteration takes
	timeStarted<-Sys.time()
	
	# Parse the specification for this set of replications

	specification<-sapply(unlist(strsplit(line, split="\t")),as.numeric)

	
	# Setting colnames will help in debugging
	names(specification)<-c(1:length(specification))
	
	# First 11 elements are scalars describing the parameters for the
	# replication. 
	
	counter<-specification[1]

	debugPrint(paste("Start of reduce task",counter))

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
	
	# Temporarily for debugging just run two models

	for(testedModelIndex in 1:testedModelCount){
		results[[testedModelIndex]]<-list()
		testedModel<-matrix(specification[(12+matrixLength*(testedModelIndex+1)):(11+matrixLength*(testedModelIndex+2))],ncol=sqrt(matrixLength))
		
		# Model estimation requires that the names are correctly included
		colnames(testedModel)<-paste("C",c(1:constructCount),sep="")
		rownames(testedModel)<-paste("C",c(1:constructCount),sep="")
		
		# Write code for testing the model

		debugPrint("Start regression")
		
		results[[testedModelIndex]]$regression <- estimateWithRegression(testedModel,data$indicators)
		
		debugPrint("Start pls")

		# TODO: Consider how well the bootstrapped construct scores correlate
		
		results[[testedModelIndex]]$plspm <- estimateWithPlspm(testedModel,data$indicators)
		
		# semPLS is very slow compared to plspm. The initial plan was to run the
		# PLS analyses twice, but to save computing time, this is now dropped
		# The line is here as a reminder, but the funtion that it calls has not
		# been really tested
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

	debugPrint("Calculate construct means")

	# There is probably a more elegant way to do this than using two nested 
	# loops.
	
	constructMeans<-list()
	for(modelTypeIndex in 1:length(results[[1]])){

		constructSums<-NULL
		
		for(testedModelIndex in 1:testedModelCount){
			
			#The results can be null if pls did not converge
			
			if(! is.null(results[[testedModelIndex]][[modelTypeIndex]])){
				if(is.null(constructSums)){
					constructSums<-results[[testedModelIndex]][[modelTypeIndex]]$constructs
				}
				else{
					constructSums<-constructSums+results[[testedModelIndex]][[modelTypeIndex]]$constructs
				}
			}
		}

		constructMeans[[modelTypeIndex]]=constructSums/testedModelCount
	}
	
	debugPrint("Start reporting")

	# Print the results. Start with the full specification string (input line)
	
	cat(line,"\n",sep="")

	#Print true factor loadings
	
	#Loop over replications and model types. Print the resutls as lines.
	
	for(modelTypeIndex in 1:length(results[[1]])){
		
		cat("\n",names(results[[1]])[modelTypeIndex],"\n")
	
		for(testedModelIndex in 1:testedModelCount){
	
			cat(names(results)[testedModelIndex],"\n")
			
			resultObj<-results[[testedModelIndex]][[modelTypeIndex]]

			#Only print things if the model converged
			
			if(is.null(resultObj)){
				cat("No convergence\n")
			}
			else{
				#Print path coefficients
				debugPrint("Paths")
				write.table(resultObj$paths,sep="\t",row.names=FALSE,col.names=FALSE)
	
	
				#Print a correlation matrix of construct scores and true scores
				debugPrint("Construct correlations")
				write.table(cor(resultObj$constructs),sep="\t",row.names=FALSE,col.names=FALSE)
	
				#Print a correlation matrix of construct scores and true scores
				debugPrint("True score correlations")
				write.table(cor(data$constructs,resultObj$constructs),sep="\t",row.names=FALSE,col.names=FALSE)
	
				#Print item-construct correlation matrix
				debugPrint("Item-construct cross-loading matrix")
				write.table(cor(data$indicators,resultObj$constructs),sep="\t",row.names=FALSE,col.names=FALSE)
	
				#Print a vector of correlations with mean of construct estimates
				debugPrint("Construct correlations with mean construct estimates")
				write.table(t(diag(cor(constructMeans[[modelTypeIndex]],resultObj$constructs))),sep="\t",row.names=FALSE,col.names=FALSE)
			}
		}
	}
	debugPrint(paste("End of reduce task",counter,"  ",timeStarted," - ",Sys.time()))

}

close(con)
