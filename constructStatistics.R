#! /usr/bin/env Rscript
#
# This file reads the combined results in data.Rdata file and calculates a
# data frame containing construct level reliability and validity statistics
#

# Needed to calculate regression from correlation matrix (mat.regress)
library(psych)

load("data.RData")

source("parameters.R")


constructStatistics=data.frame(construct=character(0),replication=numeric(0),designNumber=numeric(0),analysis=character(0),CR=numeric(0),AVE=numeric(0),minFactorLoading=numeric(0),meanFactorLoading=numeric(0),maxCrossLoading=numeric(0),maxCorrelationWithOtherConstruct=numeric(0),trueScoreCorrelation=numeric(0),deltaR2=numeric(0),sdByModels=numeric(0),sdByData=numeric(0), incomingPathsCorrect=numeric(0), incomingPathsExtra=numeric(0), incomingPathsOmitted=numeric(0), outgoingPathsCorrect=numeric(0), outgoingPathsExtra=numeric(0), outgoingPathsOmitted=numeric(0), estimatedR2=numeric(0), trueR2=numeric(0))

# Loop over all the results and calculate all needed statistics.

for(replication in 1:replications){
	for(designNumber in 1:729){
		for(analysis in 1:4){
			thisCorrelations <- correlations[[replication,designNumber,analysis]]

			# If we are dealing with partial data, do the following 
			# conditionally
			
			if(! is.null(thisCorrelations)) {
			
				constructs<-numberOfConstructs[designMatrix[designNumber,1]]
				indicators<-indicatorCounts[designMatrix[designNumber,7]]
				
				correctModel<-models[[replication,designNumber,1]]
				testedModel<-models[[replication,designNumber,2]]
				allIndicatorCols<-c((2*constructs+1):(2*constructs+constructs*indicators))
				
				for(construct in 1:constructs){
					
					print(paste("Replication:",replication,"Design number:",designNumber,"Analysis:", analysis,"Construct:",construct))
			
					# Check which row is about this construct
					
					row<-which(rownames(thisCorrelations)==paste("C",construct,sep=""))
					
					# The correlation matrix is symmetric. This assignment makes 
					# the rest of the code more readable.
				
					thisConstructCol<-row
					
					
					# Composite reliability and AVE
					#
					# Fornell, C., and Larcker, D. F. 1981. “Evaluating 
					# structural equation models with unobservable variables and 
					# measurement error,” Journal of marketing research (18:1),
					# pp. 39–50. Page 45-46
					
					
										
					indicatorCols<-which(colnames(thisCorrelations) %in% paste("i",((construct-1)*indicators+1):(construct*indicators),sep=""))

					CR<-sum(thisCorrelations[row,indicatorCols])^2/(sum(thisCorrelations[row,indicatorCols])^2+sum(1-thisCorrelations[row,indicatorCols]^2))
					
					# The denominator can be simplified since the indicators are
					# standardized
					
					AVE<-mean(thisCorrelations[row,indicatorCols]^2)
					
					# Minimum factor loading, Mean factor loading
					
					minFactorLoading<-min(thisCorrelations[row,indicatorCols])
					meanFactorLoading<-mean(thisCorrelations[row,indicatorCols])

					# Maximum cross-loading
					crossLoadingCols<-allIndicatorCols[which(!( allIndicatorCols %in% indicatorCols))]
					
					maxCrossLoading<-max(thisCorrelations[row,crossLoadingCols])
					
					# Max correlation with other construct
					
					allConstructCols<-which(colnames(thisCorrelations) %in% paste("C",1:constructs,sep=""))
					otherConstructCols<-allConstructCols[which(!( allConstructCols %in% thisConstructCol))]
					
					maxCorrelationWithOtherConstruct<-max(thisCorrelations[row,otherConstructCols])
					
					
					# Correlation with true score
					
					trueScoreCol<-which(colnames(thisCorrelations)==paste("T",construct,sep=""))
					
					trueScoreCorrelation<-thisCorrelations[row,trueScoreCol]

					# Delta R2 when other true scores added as predictors
					
					deltaR2=mat.regress(thisCorrelations,1:constructs,thisConstructCol)$R2-trueScoreCorrelation^2
					
					# Within data sd
					
					sdByData<-constructEstimateSdsByData[[replication,designNumber,analysis]][[construct]]
					
					# Within model sd
					sdByModels<-constructEstimateSdsByModel[[replication,designNumber,analysis]][[construct]]

					# Incoming paths
					
					incomingPathsCorrect <- sum(testedModel[construct,]&correctModel[construct,],na.rm = TRUE)
					incomingPathsExtra <- sum(testedModel[construct,]& ! correctModel[construct,],na.rm = TRUE)
					incomingPathsOmitted <- sum((! testedModel[construct,]) & correctModel[construct,],na.rm = TRUE)
					
					
					# Outgoing paths

					outgoingPathsCorrect <- sum(testedModel[,construct]&correctModel[,construct],na.rm = TRUE)
					outgoingPathsExtra <- sum(testedModel[,construct]& ! correctModel[,construct],na.rm = TRUE)
					outgoingPathsOmitted <- sum((! testedModel[,construct]) & correctModel[,construct],na.rm = TRUE)

					if( incomingPathsCorrect>0 | incomingPathsExtra>0 ){

						# Estimated R2 for endogenous constructs
						indices<-c(which(testedModel[construct,]==1)+constructs,constructs+construct)

						estimatedR2 <- mat.regress(thisCorrelations[indices,indices],c(1:length(indices)-1),length(indices))$R2
						
						indices<-indices-constructs
						
						# True R2 for endogenous constructs (same regression 
						# with true scores)
						
						trueR2 <- mat.regress(thisCorrelations[indices,indices],c(1:length(indices)-1),length(indices))$R2
						
					
					}
					else{
						estimatedR2<-NA
						trueR2<-NA
					}
					
					# Append the entire set to the data
					
					newRow<-data.frame(construct=paste("C",construct,sep=""),replication=replication,designNumber=designNumber,analysis=analysisTypes[analysis],CR=CR,AVE=AVE,minFactorLoading=minFactorLoading,meanFactorLoading=meanFactorLoading,maxCrossLoading=maxCrossLoading,maxCorrelationWithOtherConstruct=maxCorrelationWithOtherConstruct,trueScoreCorrelation=trueScoreCorrelation,deltaR2=deltaR2,sdByModels=sdByModels,sdByData=sdByData,incomingPathsCorrect=incomingPathsCorrect, incomingPathsExtra=incomingPathsExtra, incomingPathsOmitted=incomingPathsOmitted, outgoingPathsCorrect=outgoingPathsCorrect, outgoingPathsExtra=outgoingPathsExtra, outgoingPathsOmitted=outgoingPathsOmitted, estimatedR2=estimatedR2, trueR2=trueR2)
	
					constructStatistics<-rbind(constructStatistics,newRow)
				}
			}
		}
	}
}

print(constructStatistics)

save(constructStatistics
,file="constructStatistics.RData")



