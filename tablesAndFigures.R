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
	CR="CR",
	AVE="Root AVE",
	minFactorLoading="Minumum factor loading",
	meanFactorLoading="Mean factor loading",
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


#Only draw the tables and figures that do not yet exist.

# Figures 1, 2, and 3 are drawn with PowerPoint

######## FIGURE 4 ############

print("Figure 4")

if(!file.exists("results/figure4.pdf")){

	library(lavaan)
	
	#
	# The code below creates a set of simulated data and runs PLS, SEM, and SummedScales on these.
	#
	
	# START OF TEST PARAMETERS
	
	replications<-300
	sample<-100
	indicatorcount<-4
	factorloading<-.7
	
	# END OF TEST PARAMETERS. 
	
	#Initialize a data frame for the results
	results<-data.frame(PLS=as.numeric(NA),SEM=as.numeric(NA),Regression=as.numeric(NA),Beta=as.numeric(NA))[rep(NA,replications),]
	
	#Make a sequence from 0 to 1 consisting of uniformly distributed elements. These are the latent regression coefficients in our tests.
	
	betas<-unlist((sequence(replications)-1)/(replications-1))
	
	for(replication in 1:replications){
	
		print(paste("Running replication",replication))
	
		# GENERATE TEST DATA
	
		#Take the exogeneous latent variable from random distribution
		A<-rnorm(sample)
		
		#Calculate the values for the endogenous latent variable
	
		B<-A*betas[[replication]]+rnorm(sample)*sqrt(1-betas[[replication]]^2)
	
		#Calculate indicator values and add these to data frame with meaningful names
	
		indicators<-data.frame(row.names =c(1:sample))
		
		for ( i in 1:indicatorcount){
			indicators<-data.frame(indicators,A*factorloading+rnorm(sample)*sqrt(1-factorloading^2))
			names(indicators)[[i]]<-paste("a",i,sep="")
		}
		for ( i in 1:indicatorcount){
			indicators<-data.frame(indicators,B*factorloading+rnorm(sample)*sqrt(1-factorloading^2))
			names(indicators)[[indicatorcount+i]]<-paste("b",i,sep="")
	
		}
	
		#Run PLS 
	
		innermodel<-array(c(0,1,0,0),c(2,2))
	
		outermodel <- list(1:indicatorcount,(1:indicatorcount)+indicatorcount)
		modes = rep("A",2)
	
		tryCatch(
			pls <- plspm(indicators,innermodel,outermodel, modes)
				,correlationError = function(e){
				print(e)
				traceback()
				continue<-1
		})
		
		#Run SEM
		
		model <-paste(' B ~ A
A =~ a',paste(1:indicatorcount,collapse=" + a"),'
B =~ b',paste(1:indicatorcount,collapse=" + b"),sep="") 
		
		sem <- sem(model,data=indicators)
		
		#Calculate summed scales and run regression
	
		sumscaledata<-data.frame(A=apply(indicators[1:indicatorcount],1,mean),B=apply(indicators[(indicatorcount+1):(indicatorcount*2)],1,mean))
		
		
		regression<-lm(B~A,data=sumscaledata)
		
		#Append results to data
		results$PLS[[replication]]<-pls$path.coefs[2]
		results$SEM[[replication]]<-inspect(sem,what="std.coef")$beta[2,1]
		results$SumScale[[replication]]<-regression$coefficient[[2]]
		results$Beta[[replication]]<-betas[[replication]]
	
	}
	
	
	attach(results)
	
	pdf(file="results/figure4.pdf",height=4)
	
	par(mfrow=c(1,3)) 
	
	plot(PLS ~ Beta, asp=1, ylim=c(-.4,1.1), xlim=c(0,1))
	abline(0, 1)
	plot(SumScale ~ Beta, asp=1, ylim=c(-.4,1.1), xlim=c(0,1),ylab="Summed scales and regression")
	abline(0, 1)
	plot(SEM ~ Beta, asp=1, ylim=c(-.4,1.1), xlim=c(0,1))
	abline(0, 1)
	dev.off()
}

#Only read in this data if it is not in memory already. It takes a while to read

print("ModelData")
if( ! exists("modelData")){
	modelData <- read.delim("data/models.csv")
	print(summary(modelData))
}
print("ConstructData")

if( ! exists("constructData")){
	constructData <- read.delim("data/constructs.csv")
	
	#
	# How frequently the correlation with true score is negative
	#
	# This does not go into the table to save space. Just print to output 
	
	for(i in 1:length(analysisTypes)){
		count<-sum(constructData$analysis==i)
		negatives<-sum(constructData$analysis==i & constructData$trueScoreCorrelation<0) 
		print(paste("Analysis:",analysisTypes[i]," share of negative truescore correlations ",negatives,"/",count,"=",negatives/count))
	}
	
	constructData$trueScoreCorrelation<-abs(constructData$trueScoreCorrelation)
	constructData$minFactorLoading<-abs(constructData$minFactorLoading)
	constructData$meanFactorLoading<-abs(constructData$meanFactorLoading)
	constructData$maxCrossLoading<-abs(constructData$maxCrossLoading)

	# We typically examine the square root of AVE

	constructData$AVE<-sqrt(constructData$AVE)
	constructData$AVEMinusMaxCorrelation<-constructData$AVE-constructData$maxCorrelationWithOtherConstruct
	


}


# Table 1 and Table 2 are specified manually

# PLS is allocated the analysis types 4, 5, and 6. These are the different 
# bootstrapping options. We use only analysis types 1-4 before we discuss
# standard errors and type I and type II error rates

######## TABLE 3 ############
print("Table 3")
if(!file.exists("results/table3_full.tex")){
	# General information about quality of measurement
	
	constructDataPlus<-merge(constructData,modelData,by=c("replication","analysis","designNumber"))
	
	writeDescriptivesTable(constructDataPlus,variables=c("CR","AVE","minFactorLoading","meanFactorLoading","maxCrossLoading","AVEMinusMaxCorrelation","GlobalGoF","meanSquareResiduals","SRMR"),file="table3",analysisTypes=analysisTypes[1:4],labels=labels)
	
	rm(constructDataPlus)
}

######## TABLE 4 ############
print("Table 4")

if(!file.exists("results/table4_full.tex")){

	# Construct score reliablity and validity
	

	writeComparisonTable(constructData,variables=c("trueScoreCorrelation","deltaR2"),file="table4",analysisTypes=analysisTypes[1:4],labels=labels)
	
	# Show the experimental conditions in which PLS was best
}
print("Table 5")

######## TABLE 5 ############

if(!file.exists("results/table5_full.tex")){

	# Construct score stability
	
	writeComparisonTable(constructData,variables=c("sdByData","sdByModels"),file="table5",analysisTypes=analysisTypes[1:4],labels=labels)
}

#Only read in this data if it is not in memory already. It takes a while to read

print("RelationshipData")

if( ! exists("relationshipData")){
	relationshipData <- read.delim("data/relationships.csv")
	relationshipData$regressionARE<-abs(relationshipData$regressionTrueScore-relationshipData$regressionEstimate)

	relationshipData$correlationBias<-relationshipData$estimatedCorrelation-relationshipData$trueCorrelation*relationshipData$correlationAttenuationCoefficient
	relationshipData$correlationError<-abs(relationshipData$estimatedCorrelation-relationshipData$trueCorrelation)

	tempMatrix<-data.frame(designNumber=1:729,designMatrix)
	relationshipData<-merge(relationshipData,tempMatrix[,c("designNumber","sampleSize")])
	
	relationshipData$sampleSize<-sampleSizes[relationshipData$sampleSize]
	relationshipData$t<-relationshipData$regressionEstimate/relationshipData$regressionSE
	relationshipData$p<-(1-pt(abs(relationshipData$t),relationshipData$sampleSize-1))*2
	
	# Type I error is when the population regression coefficient is either zero 
	# or of different sign than the estimated coefficient
	
	relationshipData$TypeI05<-relationshipData$p<=.05 & relationshipData$regressionEstimate*relationshipData$regressionTrueScore <=0
	
	relationshipData$TypeI01<-relationshipData$p<=.01 & relationshipData$regressionEstimate*relationshipData$regressionTrueScore <=0

	# Type II error is when the population regression coefficient is different from zero but the p-value is non-significant sign than the estimated coefficient

	relationshipData$TypeII05<-relationshipData$p>.05 & relationshipData$regressionTrueScore == 0

	relationshipData$TypeII01<-relationshipData$p>.01 & relationshipData$regressionTrueScore == 0

	# This is not strictly speaking the right way to determine if a path is 
	# speficied because also non-convergent analyses result in NA estimates. 
	# A more approriate way would be to read the tested model from MapReduce 
	# input data or to test if the estimate is NA for all analysis types.
	

	relationshipData$specifiedAsPath<-!is.na(relationshipData$regressionEstimate)
	

	
}

######## TABLE 6 ############
print("Table 6")
if(!file.exists("results/table6_full.tex")){

	# Correlations (Hypothesis 2)
	
	# Attenuation and correlationBias

	writeComparisonTable(relationshipData,variables=c("correlationAttenuationCoefficient","correlationBias","correlationError"),file="table6",analysisTypes=analysisTypes[1:4],labels=labels)
}

#
# TODO: Check the stability of correlations across models and across data
#
print("Table 7")

######## TABLE 7 ############

if(!file.exists("results/table7_full.tex")){

	# Regression coefficients (Hypothesis 3)
	
	# Precision and accuracy
	writeComparisonTable(relationshipData,variables=c("regressionARE","regressionSE"),file="table7",analysisTypes=analysisTypes,labels=labels)
}

######## TABLE 8 ) ############

print("Table 8")


if(!file.exists("results/table8_full.tex")){

	# TYPE I and TYPE II error rate at .05
	
	writeComparisonTable(relationshipData,variables=c("TypeI05","TypeII05"),file="table8",analysisTypes=analysisTypes,labels=labels)

}

######## TABLE 9 ############

print("Table 9")

if(!file.exists("results/table9_full.tex")){

	# TYPE I and TYPE II error rate at .01

	writeComparisonTable(relationshipData,variables=c("TypeI01","TypeII01"),file="table9",analysisTypes=analysisTypes,labels=labels)

}


print("Tables 10-11")

######## TABLES 11-12 ############

#
# Two regression tables. First one has the non-bootstrapped results as 
# dependents. Second has the boostrapped results (Standard and
# Construct level changes)
#
# Independents and the design of the regressions are the same, so we do these in
# a loop.
# 

if(!file.exists("results/table10_full.tex")){

	library(lme4)

	if(!exists("dependentGroups")){
	
	# a huge regression table
	# DEPENDENTS
	# Constructs: reliability, bias
	# Correlations: attenuation, bias
	# Regressions: ARE, SE
	# Errors (.05): Type I  Type II

	# Independents
	# All experiment, construct, and relationship related covariates. Everything 
	# that has been a dependent earlier in the model
	
	designIndependents<-c("numberOfConstructs","expectedNumberOfOutgoingPaths","populationPathValues","omittedPathsShare","extraPaths","sampleSize","indicatorCount","factorLoading","factorLoadingInterval","maxErrorCorrelation","methodVariance")
	
	constructDependents<-c("trueScoreCorrelation","deltaR2")
	
	# "minFactorLoading" is dropped due to collinearity
	# TODO: Make the populationPAthValues a categorical variable.
	 constructIndependents<-c("meanFactorLoading","maxCrossLoading","trueR2","incomingPathsCorrect","incomingPathsExtra","incomingPathsOmitted","outgoingPathsCorrect","outgoingPathsExtra","outgoingPathsOmitted")
	
	correlationDependents<-c("correlationAttenuationCoefficient","correlationBias","correlationError")
	
	correlationIndependents<-c("trueCorrelation","specifiedAsPath")
	
	regressionDependents<-c("regressionARE")
	
	regressionIndependents<-c("regressionTrueScore")
	
	errorDependents<-c("regressionSE","TypeI05","TypeII05")

	relationshipIndependents=c(correlationIndependents,regressionIndependents)
	relationshipDependents=c(correlationDependents,regressionDependents)

	#
	# Create three new datasets 
	#
	constructsForRegression<-constructData[constructData$analysis==4,c("designNumber","replication","construct",constructIndependents)]

	#Set TrueR2 to zero for exogenous variables
	
	constructsForRegression[is.na(constructsForRegression[,"trueR2"]),"trueR2"]<-0

	print("Preparing construct dependents")
	
	# Take mean of all other analyses. (Each mean represents three different values)
	temp<-aggregate(constructData[constructData$analysis<4,c("designNumber","replication","construct",constructDependents)], by=list(constructData[constructData$analysis<4,"designNumber"],constructData[constructData$analysis<4,"replication"],constructData[constructData$analysis<4,"construct"]),  FUN=mean, na.rm=TRUE)
	
	temp<-merge(constructData[constructData$analysis==4,c("designNumber","replication","construct",constructDependents)],by=c("designNumber","replication","construct"),temp)
	
	temp[,constructDependents]<-temp[,paste(constructDependents,"x",sep=".")]-temp[,paste(constructDependents,"y",sep=".")]
	
	# The dependent is the difference between PLS estimate and mean of all estimates.
	constructsForRegression<-merge(constructsForRegression,temp[,c("designNumber","replication","construct",constructDependents)],by=c("designNumber","replication","construct"))

	print("Preparing relationship dependents")

	relationshipsForRegression<-relationshipData[relationshipData$analysis==4,c("designNumber","replication","to","from",relationshipIndependents)]

	temp<-aggregate(relationshipData[relationshipData$analysis<4,c("designNumber","replication","to","from",relationshipDependents)], by=list(relationshipData[relationshipData$analysis<4,"designNumber"],relationshipData[relationshipData$analysis<4,"replication"],relationshipData[relationshipData$analysis<4,"to"],relationshipData[relationshipData$analysis<4,"from"]),  FUN=mean, na.rm=TRUE)
	temp<-merge(relationshipData[relationshipData$analysis==4,c("designNumber","replication","to","from",relationshipDependents)],by=c("designNumber","replication","to","from"),temp)
	
	temp[,relationshipDependents]<-temp[,paste(relationshipDependents,"x",sep=".")]-temp[,paste(relationshipDependents,"y",sep=".")]
	
	# The dependent is the difference between PLS estimate and mean of all estimates.
	
	relationshipsForRegression<-merge(relationshipsForRegression,temp[,c("designNumber","replication","to","from",relationshipDependents)],by=c("designNumber","replication","to","from"))	


	#
	# Then do the error related things. These are a bit different since the PLS errors are estimated with three different approaches
	#

	print("Preparing error dependents")
	
	temp<-aggregate(relationshipData[relationshipData$analysis<4,c("designNumber","replication","to","from",errorDependents)], by=list(relationshipData[relationshipData$analysis<4,"designNumber"],relationshipData[relationshipData$analysis<4,"replication"],relationshipData[relationshipData$analysis<4,"to"],relationshipData[relationshipData$analysis<4,"from"]),  FUN=mean, na.rm=TRUE)


	for(i in 4:6){	
	
		temp2<-merge(relationshipData[relationshipData$analysis==i,c("designNumber","replication","to","from",errorDependents)],by=c("designNumber","replication","to","from"),temp)
	
		temp2<-cbind(temp2[,c("designNumber","replication","to","from")],temp2[,paste(errorDependents,"x",sep=".")]-temp2[,paste(errorDependents,"y",sep=".")])
		
		names(temp2)<-c("designNumber","replication","to","from",paste(errorDependents,i,sep="_"))

		temp<-merge(temp,temp2,by=c("designNumber","replication","to","from"))
	
		
	}

	relationshipsForRegression<-merge(relationshipsForRegression,temp[,c("designNumber","replication","to","from",paste(errorDependents,4,sep="_"),paste(errorDependents,5,sep="_"),paste(errorDependents,6,sep="_"))],by=c("designNumber","replication","to","from"))
	
	#
	# Merge construct related things
	#
	relationshipsForRegression<-merge(relationshipsForRegression,constructsForRegression,by.x=c("replication","designNumber","to"),by.y=c("replication","designNumber","construct"))


	relationshipsForRegression<-merge(relationshipsForRegression,constructsForRegression,by.x=c("replication","designNumber","from"),by.y=c("replication","designNumber","construct"),suffixes=c("",".from"))
	
	
	
	# Merge design related things
	
	relationshipsForRegression<-merge(relationshipsForRegression,designMatrix,by="designNumber")
	constructsForRegression<-merge(constructsForRegression,designMatrix,by="designNumber")

	print("Done merging data for regressions")
	
	# Run regressions for construct dependents
	
	# We need to adjust the cosntruct names, since C1 is not the same across replications
	
	relationshipsForRegression$to<-relationshipsForRegression$to+12*(relationshipsForRegression$replication*729+relationshipsForRegression$designNumber)
	
	relationshipsForRegression$from<-relationshipsForRegression$from+12*(relationshipsForRegression$replication*729+relationshipsForRegression$designNumber)
	
	
	
	# Recode some variables so that they make more sense in the regressions

	relationshipsForRegression[,"populationPathValues"]<-factor(relationshipsForRegression[,"populationPathValues"],labels=c("[-.5,.5]","[0,.5]","[0]"))	
	constructsForRegression[,"populationPathValues"]<-factor(constructsForRegression[,"populationPathValues"],labels=c("[-.5,.5]","[0,.5]","[0]"))

	relationshipsForRegression[,"trueCorrelation"]<-abs(relationshipsForRegression[,"trueCorrelation"])
	
	relationshipsForRegression[,"regressionTrueScore"]<-abs(relationshipsForRegression[,"regressionTrueScore"])
	
	#Standardize the data to make it easier to interpret
	
	temp<-setdiff(names(relationshipsForRegression),c("to","from","replication","desingNumber","populationPathValues","specifiedAsPath"))
	relationshipsForRegression[,temp]<-Make.Z(relationshipsForRegression[,temp])

	temp<-setdiff(names(constructsForRegression),c("replication","desingNumber","populationPathValues","specifiedAsPath"))
	constructsForRegression[,temp]<-Make.Z(constructsForRegression[,temp])

	print("Done all data preparations for regressions")
	
	

	dependentGroups<-list(constructDependents,correlationDependents,regressionDependents,c(paste(errorDependents,4,sep="_"),paste(errorDependents,5,sep="_"),paste(errorDependents,6,sep="_")))

	}
	
	print(dependentGroups)
	
	if(!exists("modelResults")){
	modelResults<-NULL

	for(dependentGroup in 1:length(dependentGroups)){
		if(dependentGroup==1){
			data<-"constructsForRegression"
		}
		else{
			data<-"relationshipsForRegression"
		}
		
		dependents<-dependentGroups[[dependentGroup]]
		
		for(dependent in 1:length(dependents)){

			# One model with experimental conditions only and one model with everything

			for(modelIndex in 1:2){
				print(paste(dependentGroup,dependent,modelIndex))
				independents <- designIndependents

				if(modelIndex>1){ 	
					independents<-c(independents,constructIndependents)
					if(dependentGroup>1) independents<-c(independents, paste(constructIndependents,"from",sep="."))
					if(dependentGroup>2) independents<-c(independents,correlationIndependents)
					if(dependentGroup>3) independents<-c(independents,regressionIndependents)
				}
				
				# All regression and errorthings have "specifiedAsPath" as always positive, so it needs to be dropped 
				if(dependentGroup>2){
					independents<-setdiff(independents,"specifiedAsPath")
				}
				
				strFormula<-paste(dependents[dependent]," ~ 1 + ",paste(independents, collapse=" + ",sep=" + ")," + (1| designNumber) + (1| replication)")

				if(dependentGroup>1 ) strFormula<-paste(strFormula," + (1|to) + (1|from)")	
				
				print(strFormula)
				reg<-lm(as.formula(gsub(" \\+ \\(.*","",strFormula)),data=get(data))
				print("Normal regression")
				print(summary(reg))
				print("Variance inflation factors from normal regression")
				print(vif(reg))
				
				modelResults<-c(modelResults,striplmer(lmer(as.formula(strFormula),data=get(data))))
				print("Stored LMER results")
			}
		}
	}
	}
	tableData<-combine.output.lmer(modelResults)

	#Print the csv table as latex. 
	formattedTableData<-tableData[,1]
	
	for(i in 1:((ncol(tableData)-2)/3)){
		formattedTableData<-cbind(formattedTableData,format(round(as.numeric(tableData[,3*i]),digits=3),nsmall=3,scientific=FALSE,trim=TRUE,na.encode=TRUE))
		ns<-abs(as.numeric(tableData[,3*i+2]))<qt(0.001,length(modelResults[[i]]@resid),lower.tail=FALSE)
		ns[is.na(ns)]<-FALSE
		formattedTableData[ns,i+1]<-paste(formattedTableData[ns,i+1],"$^{ns}$",sep="")
	}

	
	formattedTableData[formattedTableData=="NA"]<-NA
	#Remove blanck lines
	formattedTableData<-formattedTableData[!is.na(formattedTableData[1,]),]
	
	
	#Remove the intercept since it is zero
	formattedTableData<-formattedTableData[2:nrow(formattedTableData),]
	
	evenCols<-formattedTableData[1:12,1:9*2]
	formattedTableData<-formattedTableData[,1:10*2-1]
	formattedTableData[1:12,2:10]<-evenCols
	
	add.to.row<-list(list(nrow(formattedTableData)-6,nrow(formattedTableData)-7),c("\\midrule ","\\midrule "))

	# The first column contains the variable names
	helpindex=length(c(constructDependents,correlationDependents,regressionDependents))+1
	file<-"table10"
	cols<-1:helpindex
	
	print("table 1")
	print(xtable(formattedTableData[,cols]),file=paste("results/",file,"_full.tex",sep=""),add.to.row=add.to.row,sanitize.text.function=function(x){return(x)})

	print("table 2")
	print(xtable(formattedTableData[,cols]),file=paste("results/",file,"_body.tex",sep=""),include.rownames=FALSE,hline.after=NULL,only.contents=TRUE,include.colnames=FALSE,add.to.row=add.to.row,sanitize.text.function=function(x){return(x)})

	file<-"table11"
	cols<-c(1,(helpindex+2):ncol(formattedTableData))
	print(xtable(formattedTableData[,cols]),file=paste("results/",file,"_full.tex",sep=""),add.to.row=add.to.row,sanitize.text.function=function(x){return(x)})
	print(xtable(formattedTableData[,cols]),file=paste("results/",file,"_body.tex",sep=""),include.rownames=FALSE,hline.after=NULL,only.contents=TRUE,include.colnames=FALSE,add.to.row=add.to.row,sanitize.text.function=function(x){return(x)})


}


#
# Plot comparing the three PLS bootstrapping methods
#

if(!file.exists("results/figure5.pdf")){


	#
	# The code below creates a set of simulated data and runs PLS, SEM, and SummedScales on these.
	#
	
	# START OF TEST PARAMETERS
	
	replications<-50
	sample<-100
	indicatorcount<-4
	factorloading<-.7
	
	# END OF TEST PARAMETERS. 
	
	#Initialize a data frame for the results
	results<-data.frame(PLS=as.numeric(NA),sest=as.numeric(NA),seco=as.numeric(NA),sein=as.numeric(NA),Beta=as.numeric(NA))[rep(NA,replications),]
	
	#Make a sequence from 0 to 1 consisting of uniformly distributed elements. These are the latent regression coefficients in our tests.
	
	betas<-unlist((sequence(replications)-1)/(replications-1))*.5
	
	for(replication in 1:replications){
	
		print(paste("Running replication",replication))
	
		# GENERATE TEST DATA
	
		#Take the exogeneous latent variable from random distribution
		A<-rnorm(sample)
		
		#Calculate the values for the endogenous latent variable
	
		B<-A*betas[[replication]]+rnorm(sample)*sqrt(1-betas[[replication]]^2)
	
		#Calculate indicator values and add these to data frame with meaningful names
	
		indicators<-data.frame(row.names =c(1:sample))
		
		for ( i in 1:indicatorcount){
			indicators<-data.frame(indicators,A*factorloading+rnorm(sample)*sqrt(1-factorloading^2))
			names(indicators)[[i]]<-paste("a",i,sep="")
		}
		for ( i in 1:indicatorcount){
			indicators<-data.frame(indicators,B*factorloading+rnorm(sample)*sqrt(1-factorloading^2))
			names(indicators)[[indicatorcount+i]]<-paste("b",i,sep="")
	
		}
	
		#Run PLS 
	
		innermodel<-array(c(0,1,0,0),c(2,2))
	
		outermodel <- list(1:indicatorcount,(1:indicatorcount)+indicatorcount)
		modes = rep("A",2)
	
		tryCatch(
			pls <- plspm(indicators,innermodel,outermodel, modes,boot.val=TRUE)
				,correlationError = function(e){
				print(e)
				traceback()
				continue<-1
		})
		
		results$PLS[[replication]]<-pls$path.coefs[2]
		results$sest[[replication]]<-pls$boot.all$Standard$paths$Std.Error
		results$seco[[replication]]<-pls$boot.all$ConstructLevelChanges$paths$Std.Error
		results$sein[[replication]]<-pls$boot.all$IndividualSignChanges$paths$Std.Error
		results$Beta[[replication]]<-betas[[replication]]
	
	}
	
	
	attach(results)
	
	pdf(file="results/figure5.pdf",height=4)
	
	par(mar=c(5, 1, 4, 1))
	par(mfrow=c(1,3)) 
	
	errbar(Beta,PLS,PLS+2*sest,PLS-2*sest,ylab="",xlab="",yaxt="n" ,ylim=c(-0.3,0.8))
	title(main="Standard")
	abline(0, 1)

	errbar(Beta,PLS,PLS+2*seco,PLS-2*seco,ylab="",xlab="", ylim=c(-0.3,0.8),yaxt="n")
	title(main="Construct level")
	abline(0, 1)

	errbar(Beta,PLS,PLS+2*sein,PLS-2*sein,ylab="",xlab="",ylim=c(-0.3,0.8),yaxt="n")
	title(main="Indicator level")
	abline(0, 1)

	dev.off()
}

####### DISTRIBUTION PLOTS ########

#
# Show that the distribution of parameter estimates around the true value do not # follow the t-distribution. Use only data without method variance here and
# include correctly specified models. Only data with 100 observations
#

if(!file.exists("results/figure5.pdf")){
	
	# Do four of distribution plots
	
	pdf(file="results/figure5.pdf")
	
	
	par(mfrow=c(2,2)) 
	
	dataSample<-constructData
	for(i in 1:4){

		tempData <- dataSample[dataSample$analysis==i,]

		samp<-sample(1:nrow(tempData),1000)
		plot(tempData[samp,]$trueConstructR2,tempData[samp,]$estimatedConstructR2,xlab=expression(paste(R^2," for true scores")),ylab=expression(paste(R^2," for construct estimates")),main=labels[[analysisTypes[i]]])
		abline(0, 1)

	}
	dev.off()
}


#
# This figure is currently not used, so not calculated
#

if(!file.exists("results/figure6.pdf") & FALSE){
	
	# Do a distribution plot for p-value when there is no effect
	
	pdf(file="results/figure6.pdf")
	
	
	par(mfrow=c(2,2)) 
	
	#Only use data where the relationship is very close to zero
	dataSample<-relationshipData[abs(relationshipData$trueCorrelation)<.01,]
	
	tempMatrix<-data.frame(designNumber=1:729,designMatrix)
	
	useDesigns<-tempMatrix[tempMatrix$methodVariance==1,]$designNumber

	# Only use data without method variance
	
	dataSample<-dataSample[dataSample$designNumber %in% useDesigns,]
	
	# Remove NAs
	dataSample<-dataSample[!is.na(dataSample$p),]
	
	for(i in 1:length(analysisTypes)){

		tempData <- as.vector(dataSample[dataSample$analysis==i,"p"])
		tempData <- tempData[order(tempData)]
		
		x<-((1:length(tempData)))/(length(tempData))

		plot(x,tempData,xlab="Cumulative propability",ylab="Estimated value",main=labels[[analysisTypes[i]]],log="xy")
		abline(0,1)

	}
	dev.off()
}


#
# Do a similar plot as the previous except that we check if the how many 
# estimated standard errors there are between the true value and the estimate
# and what is the probablity of getting the same from drawing from the
# t-distribution.
#

#
# TODO: Make this a QQ plot
#

# There is something wrong with this figure. 

if(FALSE & !file.exists("results/figure7.pdf")){
	
	# Do a distribution plot for p-value when there is no effect
	
	#pdf(file="results/figure7.pdf")
	
	
	par(mfrow=c(2,2)) 
	
	useDesigns<-designMatrix[designMatrix[,"populationPathValues"]==3 & designMatrix[,"methodVariance"]==1 ,"designNumber"]
	
	# Only use data without method variance and the design with null effects
	dataSample<-relationshipData[relationshipData$designNumber %in% useDesigns,c("regressionEstimate","regressionTrueScore","regressionSE","analysis")]
	
	print(summary(dataSample))
	

	# Calculate how many standard errors the true score is from the estimated score
	dataSample$stat=(dataSample$regressionEstimate-dataSample$regressionTrueScore)/dataSample$regressionSE
	
	#Reverse cases where the regressionTrueScore is smaller than zero
	dataSample$stat<-dataSample$stat*abs(dataSample$regressionTrueScore)/dataSample$regressionTrueScore
	
	
	
	
	print(summary(dataSample))
	
	x <- seq(-1, 1, length=100)
	
	for(i in 1:length(analysisTypes)){

		tempData <- as.vector(dataSample[dataSample$analysis==i,"regressionEstimate"])
		
		plot(density(tempData,na.rm=TRUE),main=labels[[analysisTypes[i]]])
		lines(x, dt(x,100),col="lightgray")
	}
	#dev.off()
}

if(!file.exists("results/figure6.pdf")){
	
	# Do four of scatter plots
	
	pdf(file="results/figure6.pdf")
	
	
	par(mfrow=c(2,2)) 
	
	dataSample<-constructData
	for(i in 1:4){

		tempData <- dataSample[dataSample$analysis==i,]

		samp<-sample(1:nrow(tempData),1000)
		plot(tempData[samp,]$trueConstructR2,tempData[samp,]$estimatedConstructR2,xlab=expression(paste(R^2," for true scores")),ylab=expression(paste(R^2," for construct estimates")),main=labels[[analysisTypes[i]]])
		abline(0, 1)

	}
	dev.off()
}







