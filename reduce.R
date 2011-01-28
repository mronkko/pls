#! /usr/bin/env Rscript
# The reducer receives a population model and tested model. It runs PLS and summed scales on all specified experimental conditions for these models and outputs the results as CSV. 

source("parameters.R")

source("functions.R")

options(warn=-1)

designMatrix<-createdesignMatrix()

con <- file("stdin", open = "r")

#One line means population model- all tested models - one data specification combination. See prepare.R for details on the protocol.

while (length(line <- readLines(con, n = 1, warn = FALSE)) > 0) {

	# Record the time that this iteration takes
	timeStarted<-Sys.time()
	
	# Parse the specification for this set of replications

	specification<-sapply(unlist(strsplit(line, split="\t")),as.numeric)

	
	# First element is the counter. 
	
	counter<-specification[1]

	# Second two are the first and last row of the desing matrix 

	startIndex<-specification[2]
	endIndex<-specification[3]
	
	debugPrint(paste("Start of reduce task",counter))

	# The rest are matrices. These are population covariance, population paths,
	# and all tested models. Calculate the number and length of matrices and
	# parse the population model
	
	# Currently the tested models is always three
	testedModelCount<-3
	matrixCount<-2+testedModelCount

	matrixLength<-(length(specification)-3)/matrixCount
	constructCount<-sqrt(matrixLength)

	# Read the population model
	populationModel<-list(covariances=matrix(specification[4:(3+matrixLength)],ncol=constructCount),paths=matrix(specification[(4+matrixLength):(3+matrixLength*2)],ncol=constructCount))
	
	# Calculate construct true scores. These same scores will be used in all of 
	# the nine simulations.
	
	constructTrueScores <- populationModel$covariances,sampleSizes[sampleSizeIndex]

	# We have three data and tree models and will test all combinations. The
	# identity column for models is column 5 and the identity column for
	# indicators is column 7 of the design matrix
	
	# Test the models saving results
	
	testedModels<-list('1'=NULL,'2'=NULL,'3'=NULL)
	data<-list('1'=NULL,'2'=NULL,'3'=NULL)

	for(rowIndex in startIndex:endIndex){

		thisRow <- designMatrix[rowIndex,]
		
		# Check if model exists for this test

		if(is.null(testedModels[[thisRow[5]]])){
			testedModels[[thisRow[5]]] <- matrix(specification[(4+matrixLength*(thisRow[5]+1)):(3+matrixLength*(thisRow[5]+2))],ncol=constructCount)
		}
		
		# Check if data exists for this test 

		if(is.null(data[[thisRow[7]]])){
			data[[thisRow[7]]]<-generateData(constructTrueScores,indicatorCounts[indicatorCountIndex],factorLoadings[factorLoadingIndex],factorLoadingIntervals[factorLoadingIntervalIndex],maxErrorCorrelations[maxErrorCorrelationIndex],methodVariances[methodVarianceIndex])

		}
		
		
		# Model estimation requires that the names are correctly included
		colnames(testedModel)<-paste("C",c(1:constructCount),sep="")
		rownames(testedModel)<-paste("C",c(1:constructCount),sep="")
		
		# Write code for testing the model

		debugPrint("Start regression")
		
		results[[testedModelIndex]]$regression <- estimateWithRegression(testedModel,data$indicators)
		
		debugPrint("Start pls")

		# TODO: Consider how well the bootstrapped construct scores correlate
		
		results[[testedModelIndex]]$plspm <- estimateWithPlspm(testedModel,data$indicators)
		
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
