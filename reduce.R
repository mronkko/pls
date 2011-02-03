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

		thisDesignRow <- designMatrix[rowIndex,]

		# Check if model exists for this test. If not, read from specification.

		if(is.null(testedModels[[thisDesignRow[5]]])){
			testedModels[[thisDesignRow[5]]] <- matrix(specification[(4+matrixLength*(thisDesignRow[5]+1)):(3+matrixLength*(thisDesignRow[5]+2))],ncol=constructCount)
			# Model estimation requires that the names are correctly included
			colnames(testedModels[[thisDesignRow[5]]])<-paste("C",c(1:constructCount),sep="")
			rownames(testedModels[[thisDesignRow[5]]])<-paste("C",c(1:constructCount),sep="")

		}
		
		# Check if data exists for this test. If not, generate. 

		if(is.null(data[[thisDesignRow[7]]])){
			data[[thisDesignRow[7]]]<-generateData(constructTrueScores,indicatorCounts[thisDesignRow[7]],factorLoadings[thisDesignRow[8]],factorLoadingIntervals[thisDesignRow[9]],maxErrorCorrelations[thisDesignRow[10]],methodVariances[thisDesignRow[11]])

		}
		
		
		for(analysis in 1:length(analysisTypes)){
		
		# It is possible that the algorithm do not converge. Although this is 
		# rare, it will happen eventually. We deal with this by using try catch
		

		tryCatch({
			if(analysisTypes[analysis]=="pls"){
				results[[toString(rowIndex)]]$pls <- estimateWithPlspm(testedModels[[thisDesignRow[5]]],data[[thisDesignRow[7]]]$indicators)
			}
			else{
				results[[toString(rowIndex)]]$sumscale <- estimateWithRegression(testedModels[[thisDesignRow[5]]],data[[thisDesignRow[7]]]$indicators,method=analysisTypes[analysis])
			}},error = function(e){debugPrint(e)}
		)
			
	}
	
	#
	# Calculate estimate standard deviations
	#
	
	standardDeviations<-null
	
	for(multiplier in c(3,1)){

			for(i in 1:3){
				
				sds<-null

				# Index of the first replication that we are interested 
				# in for this calculation.
				
				baseIndex<-startIndex+(i-1)*(4-multiplier)
				index1<-toString(baseIndex)
				index2<-toString(baseIndex+multiplier)
				index3<-toString(baseIndex+2*multiplier)
				
				# Models can be con-convergent. If this happens, we do not 
				# want to calculate the sds

				if(! is.null(results[[index1]][[analysis]]) & ! is.null(results[[index2]][[analysis]]) & ! is.null(results[[index3]][[analysis]])){

					sds<-NULL
					for(constructIndex in 1:constructCount){
						
						tempMatrix<-cbind(results[[index1]][[analysis]]$constructs[,constructIndex],
						results[[index2]][[analysis]]$constructs[,constructIndex],
						results[[index3]][[analysis]]$constructs[,constructIndex])
						
						sds<-c(sds,mean(sd(t(tempMatrix))))
					}
				}
				else{
					sds<-rep(NA,constructCount)
				}
				# First three rows are within same model and the last three 
				# within data
				
				standardDeviations<-rbind(standardDeviations,sds)
			}
		}
	}
	
	#
	# After this we will analyze and report all data. 
	#
	
	# Initialize some variables that we will use in all reporting

	correctModel<- ! (is.na(populationModel$paths) | populationModel$paths==0)
	
	
	for(designNumber in startIndex:endIndex){
		for(analysis in 1:length(analysisTypes)){

			#
			# Construct level statistics first
			#
			
			thisDesingRow<-desingMatrix[desingNumber
			,]
			thisResults<-results[[toString(rowIndex)]][[analysis]]

			# If the analysis did not converge, do not write anything
			
			if(! is.null(thisResults)){
			thisCorrelations <-cor(cbind(constructTrueScores,thisResults$constructs,data[[thisDesignRow[7]]]$indicators)

			indicatorCount<-indicatorCounts[thisDesignRow[7]]

			testedModel<-models[thisDesignRow[5]]]

			# Which columns of thisCorrelations are indicator correlations
			allIndicatorCols<-c((2*constructCount+1):(2*constructCount+constructCount*indicatorCount))

			allConstructCols<-(constructCount+1):(constructCount*2)
			
			for(construct in 1:constructCount){
				
				row=constructCount+construct
				
				# The correlation matrix is symmetric. This assignment makes 
				# the rest of the code more readable.
			
				thisConstructCol<-row
				
				
				# Composite reliability and AVE
				#
				# Fornell, C., and Larcker, D. F. 1981. “Evaluating 
				# structural equation models with unobservable variables and 
				# measurement error,” Journal of marketing research (18:1),
				# pp. 39–50. Page 45-46
				
				
				temp<-2*constructCount+(construct-1)*indicatorCount
				indicatorCols<-(temp+1):(temp+indicatorCount)
				CR<-sum(thisCorrelations[row,indicatorCols])^2/(sum(thisCorrelations[row,indicatorCols])^2+sum(1-thisCorrelations[row,indicatorCols]^2))
				
				# The denominator can be simplified since the indicatorCount are
				# standardized
				
				AVE<-mean(thisCorrelations[row,indicatorCols]^2)
				
				# Minimum factor loading, Mean factor loading
				
				minFactorLoading<-min(thisCorrelations[row,indicatorCols])
				meanFactorLoading<-mean(thisCorrelations[row,indicatorCols])

				# Maximum cross-loading
				crossLoadingCols<-allIndicatorCols[which(!( allIndicatorCols %in% indicatorCols))]
				
				maxCrossLoading<-max(thisCorrelations[row,crossLoadingCols])
				
				# Max correlation with other construct
				
				otherConstructCols<-allConstructCols[which(!( allConstructCols %in% thisConstructCol))]
				
				maxCorrelationWithOtherConstruct<-max(thisCorrelations[row,otherConstructCols])
				
				
				# Correlation with true score
				
				trueScoreCol<-construct
				
				trueScoreCorrelation<-thisCorrelations[row,trueScoreCol]

				# Delta R2 when other true scores added as predictors
				
				deltaR2=mat.regress(thisCorrelations,1:constructCount,thisConstructCol)$R2-trueScoreCorrelation^2
				
				# Within data sd
				
				sdByData<-standardDeviations[thisDesingRow[5]+3,construct]
				
				# Within model sd
				
				sdByModels<-standardDeviations[thisDesingRow[7],construct]
				
				# Incoming paths
				
				incomingPathsCorrect <- sum(testedModel[construct,]&correctModel[construct,],na.rm = TRUE)
				incomingPathsExtra <- sum(testedModel[construct,]& ! correctModel[construct,],na.rm = TRUE)
				incomingPathsOmitted <- sum((! testedModel[construct,]) & correctModel[construct,],na.rm = TRUE)
				
				
				# Outgoing paths

				outgoingPathsCorrect <- sum(testedModel[,construct]&correctModel[,construct],na.rm = TRUE)
				outgoingPathsExtra <- sum(testedModel[,construct]& ! correctModel[,construct],na.rm = TRUE)
				outgoingPathsOmitted <- sum((! testedModel[,construct]) & correctModel[,construct],na.rm = TRUE)

				if( incomingPathsCorrect>0 | incomingPathsExtra>0 ){

					# Estimated R2 for endogenous constructCount
					indices<-c(which(testedModel[construct,]==1)+constructCount,constructCount+construct)

					estimatedR2 <- mat.regress(thisCorrelations[indices,indices],c(1:length(indices)-1),length(indices))$R2
					
					indices<-indices-constructCount
					
					# True R2 for endogenous constructCount (same regression 
					# with true scores)
					
					trueR2 <- mat.regress(thisCorrelations[indices,indices],c(1:length(indices)-1),length(indices))$R2
					
				
				}
				else{
					estimatedR2<-NA
					trueR2<-NA
				}
				
				# Start a new row by printing "C" to mark this as a row that is # about constructs,the replication, desingNumber,
				# constructNumber and then the results
				
				cat("C",replication,desingNumber,analysis,construct,CR, AVE, minFactorLoading, meanFactorLoading, maxCrossLoading, maxCorrelationWithOtherConstruct, trueScoreCorrelation, deltaR2, estimatedR2, trueR2, sdByData, sdByModels, incomingPathsCorrect, incomingPathsExtra, incomingPathsOmitted, outgoingPathsCorrect, outgoingPathsExtra, outgoingPathsOmitted sep="\t")
				cat("\n")
			}
	
			#
			# Then relationship level statistics
			#
	
			thisPaths<-data[[thisDesignRow[7]]]$paths
		
			# Loop over all correlations in the lower diagonal
			
			for(from in 1:(constructCount-1)){
				for(to in (from+1):constructCount){
			

					# CORRELATIONS (hypothesis 2)
					
					# Correlation is the sum of attenuation and bias. 
					# Calculate the attenuation and then we know the bias 
					# too

					trueCorrelation<-thisCorrelations[from,to]
					estimatedCorrelation<-thisCorrelations[from+constructCount,to+constructCount]
					
					# Reliability of the construct estimate is the square of 
					# the correlation between the construct and the true
					# score. Hence the attenuation coefficient is the 
					# product of the true score correlations. See Cohen p. 
					# 55-56
					
					attenuationCoefficient<-thisCorrelations[from,from]*thisCorrelations[to,to]
					
					# Bias is estimated correlation minus true correlation 
					# times attenuation. - This is a trivial calculation, so # we do not store the result
					
					# bias<-estimatedCorrelation-trueCorrelation*attenuationCoefficient
					
					
					# REGRESSION COEFFICIENTS (hypothesis 3)
					
					# Precision is the SD of difference between true 
					# regression coefficient and the estimate
					# Accuracy is the Mean difference between true 
					# regression coefficient and the estimate
					
					
			
					regressionTrueScore<-populationModel$paths[to,from]
					regressionEstimate<-thisPaths[thisPaths$To==to & thisPaths$From==from,"Estimate"]
					# STANDARD ERRORS (hypothesis 4)
					
					regressionSE<-thisPaths[thisPaths$To==to & thisPaths$From==from,"Estimate"]

					cat("R",replication,desingNumber,analysis,to,from,trueCorrelation, estimatedCorrelation, attenuationCoefficient, regressionTrueScore, regressionEstimate, regressionSE,sep="\t")
					
					cat("\n")

				}
			}
		}
	}
	debugPrint(paste("End of reduce task",counter,"  ",timeStarted," - ",Sys.time()))
}

close(con)
