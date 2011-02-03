#!/usr/bin/env Rscript
# The reducer receives a population model and tested model. It runs PLS and summed scales on all specified experimental conditions for these models and outputs the results as CSV. 

source("include/parameters.R")
source("include/functions.R")

options(warn=-1)

designMatrix<-createdesignMatrix()

con <- file("stdin", open = "r")

#One line means population model- all tested models - one data specification combination. See prepare.R for details on the protocol.

while (length(line <- trim(readLines(con, n = 1, warn = FALSE))) > 0) {

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
		
		
		debugPrint("Start regression")
		
		
		# It is possible that the algorithm do not converge. Although this is 
		# rare, it will happen eventually. We deal with this by using try catch
		

		tryCatch(
			results[[toString(rowIndex)]]$sumscale <- estimateWithRegression(testedModels[[thisRow[5]]],data[[thisRow[7]]]$indicators,method="sumscale")
			,error = function(e){debugPrint(e)}
		)

		tryCatch(
			results[[toString(rowIndex)]]$component <- estimateWithRegression(testedModels[[thisRow[5]]],data[[thisRow[7]]]$indicators,method="component")
			,error = function(e){debugPrint(e)}
		)

		tryCatch(
			results[[toString(rowIndex)]]$factor <- estimateWithRegression(testedModels[[thisRow[5]]],data[[thisRow[7]]]$indicators,method="factor")
			,error = function(e){debugPrint(e)}
		)

		
		debugPrint("Start pls")
		
		tryCatch(
			results[[toString(rowIndex)]]$pls <- estimateWithPlspm(testedModels[[thisRow[5]]],data[[thisRow[7]]]$indicators)
			,error = function(e){debugPrint(e)}
		)
			
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

	cat("\nDesign\n")
	cat(line,"\n",sep="")

	
	#Loop over the three indicator data
	
	for(i in 1:3){
		cat("\nData ",i,"\n",sep="")
		#Print true factor loadings
		cat("\nPopulation factor loadings\n")
		tempMatrix<-matrix(data[[i]]$factorLoadings,nrow=1)
		colnames(tempMatrix)<-names(data[[i]]$indicators)
		write.table(tempMatrix,sep="\t")
		
	}

	#Loop ove model types and print results

	for(modelTypeIndex in 1:length(results[[toString(startIndex)]])){
		
		cat("\n\n",names(results[[toString(startIndex)]])[modelTypeIndex],"\n\n")
	
		# Print calculate standard deviations of construct within cases over
		# each data and each model. Calculate mean of sd by construct and print 
		# to results

		cat("\nConstruct score sds over data/model\n")
		for(multiplier in c(3,1)){
			labels <- c("Model","","Data")
			label<-labels[multiplier]
			for(i in 1:3){
				# Index of the first replication that we are interested 
				# in for this calculation.
				
				baseIndex<-startIndex+(i-1)*(4-multiplier)
				index1<-toString(baseIndex)
				index2<-toString(baseIndex+multiplier)
				index3<-toString(baseIndex+2*multiplier)
				
				# PLS models can be con-convergent. If this happens, we do not 
				# want to calculate the sds

				# These are for debugging only. Need to be removed later.
				cat("\n")
				debugPrint(results[[index1]][[modelTypeIndex]])
				debugPrint(results[[index2]][[modelTypeIndex]])
				debugPrint(results[[index3]][[modelTypeIndex]])

				if(! is.null(results[[index1]][[modelTypeIndex]]) & ! is.null(results[[index2]][[modelTypeIndex]]) & ! is.null(results[[index3]][[modelTypeIndex]])){

					cat(paste("\n",label,i,"\n"))

					sds<-NULL
					for(constructIndex in 1:constructCount){
						
						tempMatrix<-cbind(results[[index1]][[modelTypeIndex]]$constructs[,constructIndex],
						results[[index2]][[modelTypeIndex]]$constructs[,constructIndex],
						results[[index3]][[modelTypeIndex]]$constructs[,constructIndex])
						
						sds<-c(sds,mean(sd(t(tempMatrix))))
					}
					tempMatrix<-matrix(sds,nrow=1)
					colnames(tempMatrix)<-names(data[[i]]$constructs)
					write.table(tempMatrix,sep="\t")
					cat("\n")
				}
			}
		}
	
		# Loop over the nine tested models and print results
		for(rowIndex in startIndex:endIndex){

			thisRow <- designMatrix[rowIndex,]

			cat("\nRow index",rowIndex)
			cat("\n\n")
			
			resultObj<-results[[toString(rowIndex)]][[modelTypeIndex]]

			#Only print things if the model converged
			
			if(is.null(resultObj)){
				cat("\nNo convergence\n")
			}
			else{
				#Print path coefficients
				cat("\nPaths\n")
				write.table(resultObj$paths,sep="\t")
	
	
				#Print a correlation matrix of construct scores, true scores, andd indicators
				cat("\nCorrelations\n")
				write.table(cor(cbind(data[[thisRow[7]]]$constructs,resultObj$constructs,data[[thisRow[7]]]$indicators)),sep="\t")
	
			}
		}
	}
	debugPrint(paste("End of reduce task",counter,"  ",timeStarted," - ",Sys.time()))

}

close(con)
