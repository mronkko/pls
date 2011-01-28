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
	
	debugPrint("Genarate construct true scores")
	
	constructTrueScores <- mvrnorm(n=sampleSizes[designMatrix[[startIndex,4]]],rep(0,constructCount),populationModel$covariances)


	# We have three data and tree models and will test all combinations. The
	# identity column for models is column 5 and the identity column for
	# indicators is column 7 of the design matrix
	
	# Test the models saving results
	
	testedModels<-list('1'=NULL,'2'=NULL,'3'=NULL)
	data<-list('1'=NULL,'2'=NULL,'3'=NULL)
	results<-list()
	
	for(rowIndex in startIndex:endIndex){

		thisRow <- designMatrix[rowIndex,]

#		debugPrint(paste("Row ",rowIndex,"of desing matrix"))
#		debugPrint(thisRow)
		
		# Check if model exists for this test. If not, read from specification.

		if(is.null(testedModels[[thisRow[5]]])){
			testedModels[[thisRow[5]]] <- matrix(specification[(4+matrixLength*(thisRow[5]+1)):(3+matrixLength*(thisRow[5]+2))],ncol=constructCount)
			# Model estimation requires that the names are correctly included
			colnames(testedModels[[thisRow[5]]])<-paste("C",c(1:constructCount),sep="")
			rownames(testedModels[[thisRow[5]]])<-paste("C",c(1:constructCount),sep="")

		}
		
		# Check if data exists for this test. If not, generate. 

		if(is.null(data[[thisRow[7]]])){
			data[[thisRow[7]]]<-generateData(constructTrueScores,indicatorCounts[thisRow[7]],factorLoadings[thisRow[8]],factorLoadingIntervals[thisRow[9]],maxErrorCorrelations[thisRow[10]],methodVariances[thisRow[11]])

		}
		
#		debugPrint(testedModels[[thisRow[5]]])
#		debugPrint(summary(data[[thisRow[7]]]$indicators))
		
		debugPrint("Start regression")
		
		results[[toString(rowIndex)]]$sumscale <- estimateWithRegression(testedModels[[thisRow[5]]],data[[thisRow[7]]]$indicators,method="sumscale")

		results[[toString(rowIndex)]]$component <- estimateWithRegression(testedModels[[thisRow[5]]],data[[thisRow[7]]]$indicators,method="component")

		results[[toString(rowIndex)]]$factor <- estimateWithRegression(testedModels[[thisRow[5]]],data[[thisRow[7]]]$indicators,method="factor")

		
		debugPrint("Start pls")

		results[[toString(rowIndex)]]$pls <- estimateWithPlspm(testedModels[[thisRow[5]]],data[[thisRow[7]]]$indicators)
		
	}
	
	# Calculate some statistics that need to be calculated over all analyses
	# done for this particular set of replications
	
	#
	# Calculate correlations for each construct over all estimated models.
	# We do this by calculating a mean score over all nine tests for each type
	# of test and then test calculating how well each replication correlates 
	# with the mean values
	#


	debugPrint("Start reporting")

	# Print the results. Start with the full specification string (input line)

	cat("Design\n")
	cat(line,"\n",sep="")

	cat("Data\n")
	#Print construct sample correlations
	write.table(cor(constructTrueScores),sep="\t",row.names=FALSE,col.names=FALSE)
	
	#Loop over the three indicator data
	
	for(i in 1:3){

		#Print true factor loadings
		cat(data[[i]]$factorLoadings,sep="\t")
		cat("\n")

		#Print correlation matrix of items
		write.table(cor(data[[i]]$indicators),sep="\t",row.names=FALSE,col.names=FALSE)

		#Print item-construct correlations
		write.table(cor(data[[i]]$indicators,constructTrueScores),sep="\t",row.names=FALSE,col.names=FALSE)
		
	}

	#Loop ove model types and print results

	for(modelTypeIndex in 1:length(results[[toString(startIndex)]])){
		
		cat("\n",names(results[[toString(startIndex)]])[modelTypeIndex],"\n")
	
		# Print calculate standard deviations of construct within cases over
		# each data. Calculate mean of sd by construct and print to results

		for(dataIndex in 1:3){
			for(constructIndex in 1:constructCount){
				
				tempMatrix<-cbind(results[[toString(startIndex+dataIndex-1)]][modelTypeIndex]$constructs[,constructIndex],
				results[[toString(startIndex+dataIndex+3)]][[modelTypeIndex]]$constructs[,constructIndex],
				results[[toString(startIndex+dataIndex+6)]][[modelTypeIndex]]$constructs[,constructIndex])
				
				if(constructIndex!=1) cat("\t")
				cat(mean(sd(t(tempMatrix))))
			}
			cat("\n")
		}
		
		# Print calculate standard deviations of construct within cases over
		# each model. Calculate mean of sd  by construct and print to results

		
		for(modelIndex in 1:3){
			for(constructIndex in 1:constructCount){
				tempMatrix<-cbind(results[[toString(startIndex+modelIndex*3-3)]][modelTypeIndex]$constructs[,constructIndex],
				results[[toString(startIndex+modelIndex*3-2)]][[modelTypeIndex]]$constructs[,constructIndex],
				results[[toString(startIndex+modelIndex*3-1)]][[modelTypeIndex]]$constructs[,constructIndex])

				if(constructIndex!=1) cat("\t")
				cat(mean(sd(t(tempMatrix))))
			}
			cat("\n")
		}
	
		# Loop over the nine tested models and print results
		for(rowIndex in startIndex:endIndex){

			thisRow <- designMatrix[rowIndex,]

			cat(rowIndex)
			cat("\n")
			
			resultObj<-results[[toString(rowIndex)]][[modelTypeIndex]]

			#Only print things if the model converged
			
			if(is.null(resultObj)){
				cat("No convergence\n")
			}
			else{
				#Print path coefficients
				cat("Paths\n")
				write.table(resultObj$paths,sep="\t",row.names=FALSE,col.names=FALSE)
	
	
				#Print a correlation matrix of construct scores and true scores
				cat("Construct correlations\n")
				write.table(cor(resultObj$constructs),sep="\t",row.names=FALSE,col.names=FALSE)
	
				#Print a correlation matrix of construct scores and true scores
				cat("True score correlations\n")
				
				write.table(cor(data[[thisRow[7]]]$constructs,resultObj$constructs),sep="\t",row.names=FALSE,col.names=FALSE)
	
				#Print item-construct correlation matrix
				cat("Item-construct cross-loading matrix\n")
				write.table(cor(data[[thisRow[7]]]$indicators,resultObj$constructs),sep="\t",row.names=FALSE,col.names=FALSE)
	
			}
		}
	}
	debugPrint(paste("End of reduce task",counter,"  ",timeStarted," - ",Sys.time()))

}
wc
close(con)
