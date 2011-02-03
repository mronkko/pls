#! /usr/bin/env Rscript
#
# This file reads the combined results in data.Rdata file and calculates a
# data frame containing correlation and path related statistics
#

# Needed to calculate regression from correlation matrix
library(psych)

load("data.RData")

source("include/parameters.R")


relationshipStatistics=data.frame(from=character(0),to=character(0),replication=numeric(0),designNumber=numeric(0),analysis=character(0),CR=numeric(0),AVE=numeric(0),minFactorLoading=numeric(0),meanFactorLoading=numeric(0),maxCrossLoading=numeric(0),maxCorrelationWithOtherConstruct=numeric(0),trueScoreCorrelation=numeric(0),deltaR2=numeric(0),sdByModels=numeric(0),sdByData=numeric(0))

# Loop over all the results and calculate all needed statistics.

for(replication in 1:replications){
	for(designNumber in 1:729){
		thisPaths<-paths[paths$replication==replication & paths$designNumber==designNumber & paths$analysis=="truevalue",]

		for(analysis in 1:4){

			thisCorrelations <- correlations[[replication,designNumber,analysis]]

			thisPaths<-paths[paths$replication==replication & paths$designNumber==designNumber & paths$analysis==analysis,]
			
			# If we are dealing with partial data, do the following 
			# conditionally
			
			if(! is.null(thisCorrelations)) {
			
				constructs<-numberOfConstructs[designMatrix[designNumber,1]]
				indicators<-indicatorCounts[designMatrix[designNumber,7]]
				
				allIndicatorCols<-c((2*constructs+1):(2*constructs+constructs*indicators))
				
				# Loop over all correlations in the lower diagonal
				
				for(from in 1:(constructs-1)){
					for(to in (from+1):constructs){
				
						print(paste("Replication:",replication,"Design number:",designNumber,"Analysis:", analysis,"From:",from,"To:",to))

						# CORRELATIONS (hypothesis 2)
						
						# Correlation is the sum of attenuation and bias. 
						# Calculate the attenuation and then we know the bias 
						# too

						trueCorrelation<-thisCorrelations[from,to]
						estimatedCorrelation<-thisCorrelations[from+constructs,to+constructs]
						
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
						
						
						# colnames(newRow)<-c("From","To","Estimate","Mean.Boot","Std.Error","perc.05","perc.95","replication","designNumber","analysis")
				
						regressionTrueScore<-thisTruePaths[thisPaths$To==to & thisPaths$From==from,"Estimate"]
						regressionEstimate<-thisPaths[thisPaths$To==to & thisPaths$From==from,"Estimate"]
						# STANDARD ERRORS (hypothesis 4)
						
						regressionSE<-thisPaths[thisPaths$To==to & thisPaths$From==from,"Estimate"]

						tValue<-regressionEstimate/regressionSE
						pValue<-pt(abs(tValue),designMatrix[designNumber,"sampleSize"],lower.tail=FALSE)
						
	
						# Append the entire set to the data
						
						newRow<-data.frame(construct=paste("C",construct,sep=""),replication=replication,designNumber=designNumber,analysis=analysisTypes[analysis],CR=CR,AVE=AVE,minFactorLoading=minFactorLoading,meanFactorLoading=meanFactorLoading,maxCrossLoading=maxCrossLoading,maxCorrelationWithOtherConstruct=maxCorrelationWithOtherConstruct,trueScoreCorrelation=trueScoreCorrelation,deltaR2=deltaR2,sdByModels=sdByModels,sdByData=sdByData)
						constructStatistics<-rbind(constructStatistics,newRow)
					}
				}
			}
		}
	}
}

print(constructStatistics)

save(constructStatistics
,file="relationshipStatistics.RData")



