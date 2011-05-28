#! /usr/bin/env Rscript
#
# This file creates a set of simulated data and runs PLS, SEM, and SummedScales on these.
#


library(car)
library(QuantPsyc)

# Read simulation parameters
source("include/parameters.R")

# Read function definitions
source("include/functions.R")
source("include/functionsTablesAndFigures.R")


# Define labels for variables

labels<-list(
	CR="Composite reliability",
	AVE="Root AVE",
	minFactorLoading="Minumum composite loading",
	meanFactorLoading="Mean composite loading",
	maxCrossLoading="Maximum cross loading",
	AVEMinusMaxCorrelation="Root AVE - max correlation",
	sumscale="Summed scales",
	component="Components",
	factor="Factors",
	pls_Standard="PLS",
	pls_IndividualSignChanges="PLS, Indicator sign correction",
	pls_ConstructLevelChanges="PLS, Construct sign correction",
	numberOfConstructs="Number of constructs",
	expectedNumberOfOutgoingPaths="Expected number of paths",
	populationPathValues="Population path values",
	omittedPathsShare="Omitted paths",
	extraPaths="Extra paths",
	sampleSize="Sample size",
	indicatorCount="Indicators",
	factorLoading="Meand factor loading",
	factorLoadingInterval="Factor loading variation",
	maxErrorCorrelation="Max error correlations",
	methodVariance="Method variances",GlobalGoF="Global goodness of fit",meanSquareResiduals="Indicator mean square residual",SRMR="SRMR")

designMatrix<-createdesignMatrix()



#Only read in this data if it is not in memory already. It takes a while to read

print("ModelData")
if( ! exists("modelData")){
	modelData <- read.delim("data/models.csv")
	print(summary(modelData))
}
print("constructDataAppendix")

if( ! exists("constructDataAppendix")){
	constructDataAppendix <- read.delim("data/constructs.csv")
	
	
	
	#
	# How frequently the correlation with true score is negative
	#
	# This does not go into the table to save space. Just print to output 
	
	for(i in 1:length(analysisTypes)){
		count<-sum(constructDataAppendix$analysis==i)
		negatives<-sum(constructDataAppendix$analysis==i & constructDataAppendix$trueScoreCorrelation<0) 
		print(paste("Analysis:",analysisTypes[i]," share of negative truescore correlations ",negatives,"/",count,"=",negatives/count))
	}
	
	constructDataAppendix$trueScoreCorrelation<-abs(constructDataAppendix$trueScoreCorrelation)
	constructDataAppendix$minFactorLoading<-abs(constructDataAppendix$minFactorLoading)
	constructDataAppendix$meanFactorLoading<-abs(constructDataAppendix$meanFactorLoading)
	constructDataAppendix$maxCrossLoading<-abs(constructDataAppendix$maxCrossLoading)

	# We typically examine the square root of AVE

	constructDataAppendix$AVE<-sqrt(constructDataAppendix$AVE)
	constructDataAppendix$AVEMinusMaxCorrelation<-constructDataAppendix$AVE-constructDataAppendix$maxCorrelationWithOtherConstruct
	

	#Limit to the following criteria:
	#
	#Composite reliability > .7, Minumum composite loading > .7, Root AVE > .5, #Maximum cross loading < .4, and Root AVE - max correlation > 0.
	#
	
	indices<-(constructDataAppendix$minFactorLoading < .7 | constructDataAppendix$CR  < .7 | constructDataAppendix$AVE <.5 | constructDataAppendix$maxCrossLoading >.4 | constructDataAppendix$AVEMinusMaxCorrelation < 0) 
	
	invalidModels<-unique(constructDataAppendix$replication[indices]*1000+constructDataAppendix$designNumber[indices])

	constructDataAppendix<-constructDataAppendix[! ( (constructDataAppendix$replication*1000 + constructDataAppendix$designNumber) %in% invalidModels),]
}



# PLS is allocated the analysis types 4, 5, and 6. These are the different 
# bootstrapping options. We use only analysis types 1-4 before we discuss
# standard errors and type I and type II error rates

######## TABLE 3 ############
print("Table 3")
if(!file.exists("results/table3appendix_full.tex")){
	# General information about quality of measurement
	
	constructDataAppendixPlus<-merge(constructDataAppendix,modelData,by=c("replication","analysis","designNumber"))
	
	#TODO: Make meanSquareResidual needs to be taken a square root of
	
	constructDataAppendixPlus$SRMR<-sqrt(constructDataAppendixPlus$SRMR)
	writeDescriptivesTable(constructDataAppendixPlus,variables=c("CR","AVE","minFactorLoading","meanFactorLoading","maxCrossLoading","AVEMinusMaxCorrelation","meanSquareResiduals"),file="table3appendix",analysisTypes=analysisTypes[1:4],labels=labels)
	
	rm(constructDataAppendixPlus)
}

######## TABLE 4 ############
print("Table 4")

if(!file.exists("results/table4appendix_full.tex")){

	# Construct score reliablity and validity
	

	writeComparisonTable(constructDataAppendix,variables=c("trueScoreCorrelation","deltaR2"),file="table4appendix",analysisTypes=analysisTypes[1:4],labels=labels)
	
	# Show the experimental conditions in which PLS was best
}
print("Table 5")

######## TABLE 5 ############

if(!file.exists("results/table5appendix_full.tex")){

	# Construct score stability
	
	writeComparisonTable(constructDataAppendix,variables=c("sdByData","sdByModels"),file="table5appendix",analysisTypes=analysisTypes[1:4],labels=labels)
}

#Only read in this data if it is not in memory already. It takes a while to read

print("RelationshipDataAppendix")

if( ! exists("RelationshipDataAppendix")){
	RelationshipDataAppendix <- read.delim("data/relationships.csv")
	RelationshipDataAppendix$regressionARE<-abs(RelationshipDataAppendix$regressionTrueScore-RelationshipDataAppendix$regressionEstimate)

	RelationshipDataAppendix$correlationBias<-RelationshipDataAppendix$estimatedCorrelation-RelationshipDataAppendix$trueCorrelation*RelationshipDataAppendix$correlationAttenuationCoefficient
	RelationshipDataAppendix$correlationError<-abs(RelationshipDataAppendix$estimatedCorrelation-RelationshipDataAppendix$trueCorrelation)

	tempMatrix<-data.frame(designNumber=1:729,designMatrix)
	RelationshipDataAppendix<-merge(RelationshipDataAppendix,tempMatrix[,c("designNumber","sampleSize")])
	
	RelationshipDataAppendix$sampleSize<-sampleSizes[RelationshipDataAppendix$sampleSize]
	RelationshipDataAppendix$t<-RelationshipDataAppendix$regressionEstimate/RelationshipDataAppendix$regressionSE
	RelationshipDataAppendix$p<-(1-pt(abs(RelationshipDataAppendix$t),RelationshipDataAppendix$sampleSize-1))*2
	
	# Type I error is when the population regression coefficient is either zero 
	# or of different sign than the estimated coefficient
	
	RelationshipDataAppendix$TypeI05<-RelationshipDataAppendix$p<=.05 & RelationshipDataAppendix$regressionEstimate*RelationshipDataAppendix$regressionTrueScore <=0
	
	RelationshipDataAppendix$TypeI01<-RelationshipDataAppendix$p<=.01 & RelationshipDataAppendix$regressionEstimate*RelationshipDataAppendix$regressionTrueScore <=0

	# Type II error is when the population regression coefficient is different from zero but the p-value is non-significant sign than the estimated coefficient

	RelationshipDataAppendix$TypeII05<-RelationshipDataAppendix$p>.05 & RelationshipDataAppendix$regressionTrueScore == 0

	RelationshipDataAppendix$TypeII01<-RelationshipDataAppendix$p>.01 & RelationshipDataAppendix$regressionTrueScore == 0

	# This is not strictly speaking the right way to determine if a path is 
	# speficied because also non-convergent analyses result in NA estimates. 
	# A more approriate way would be to read the tested model from MapReduce 
	# input data or to test if the estimate is NA for all analysis types.
	

	RelationshipDataAppendix$specifiedAsPath<-!is.na(RelationshipDataAppendix$regressionEstimate)
	
	# Eliminate replications that do not meet the quality criteria

	RelationshipDataAppendix<-RelationshipDataAppendix[! ( (RelationshipDataAppendix$replication*1000 + RelationshipDataAppendix$designNumber) %in% invalidModels),]

	
}

######## TABLE 6 ############
print("Table 6")
if(!file.exists("results/table6appendix_full.tex")){

	# Correlations (Hypothesis 2)
	
	# Attenuation and correlationBias

	writeComparisonTable(RelationshipDataAppendix,variables=c("correlationAttenuationCoefficient","correlationBias","correlationError"),file="table6appendix",analysisTypes=analysisTypes[1:4],labels=labels)
}

#
# TODO: Check the stability of correlations across models and across data
#
print("Table 7")

######## TABLE 7 ############

if(!file.exists("results/table7appendix_full.tex")){

	# Regression coefficients (Hypothesis 3)
	
	# Precision and accuracy
	writeComparisonTable(RelationshipDataAppendix,variables=c("regressionARE","regressionSE"),file="table7appendix",analysisTypes=analysisTypes,labels=labels)
}

######## TABLE 8 ) ############

print("Table 8")


if(!file.exists("results/table8appendix_full.tex")){

	# TYPE I and TYPE II error rate at .05
	
	writeComparisonTable(RelationshipDataAppendix,variables=c("TypeI05","TypeII05"),file="table8appendix",analysisTypes=analysisTypes,labels=labels)

}

######## TABLE 9 ############

print("Table 9")

if(!file.exists("results/table9appendix_full.tex")){

	# TYPE I and TYPE II error rate at .01

	writeComparisonTable(RelationshipDataAppendix,variables=c("TypeI01","TypeII01"),file="table9appendix",analysisTypes=analysisTypes,labels=labels)

}






