#! /usr/bin/env Rscript
#
# This file creates a set of simulated data and runs PLS, SEM, and SummedScales on these.
#

print("Start of file")

#setwd("/Users/mronkko/Documents/Research/pls")

#library(car)
library(QuantPsyc)
library(Hmisc)
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
	factorLoading="Mean factor loading",
	factorLoadingInterval="Factor loading variation",
	maxErrorCorrelation="Max error correlations",
	methodVariance="Method variances",GlobalGoF="Global goodness of fit",meanSquareResiduals="Indicator mean square residual",
	SRMR="SRMR",
	trueScoreR2="Variance explained by true score",
	deltaR2constructs="Add. variance explained by other constructs",
	deltaR2errors="Add. variance explained by error terms",
	totalR2="Total variance explained",
	sdByData="SD of construct scores over three data",
	sdByModels="SD of construct scores over three models",
	regressionARE="Absolute relative error",
	TypeI05 = "False positive, p<.05",
	TypeI01 = "False positive, p<.01",
	TypeI001 = "False positive, p<.001",
	TypeII05 = "False negative, p<.05",
	TypeII01 = "False negative, p<.01",
	TypeII001 = "False negative, p<.001",
	TypeIII05 = "Sig. eff. in the wrong dir., p<.05",
	TypeIII01 = "Sig. eff. in the wrong dir., p<.01",
	TypeIII001 = "Sig. eff. in the wrong dir., p<.001"
	)

designMatrix<-createdesignMatrix()


#Only draw the tables and figures that do not yet exist.

# Figure 1 is drawn with PowerPoint


print("Figure 2")

if(!file.exists("results/figure2.pdf")){

	#
	# The code below creates a set of simulated data and runs PLS, Principal components, and SummedScales on these.
	#
	
	# START OF TEST PARAMETERS
	
	replications<-300
	sample<-100
	indicatorcount<-3
	factorloadings<-c(.7,.9)
	
	# END OF TEST PARAMETERS. 
	
	#Initialize a data frame for the results
	results<-data.frame(PLS=as.numeric(NA),SumScale=as.numeric(NA),PrinComp=as.numeric(NA),Factor=as.numeric(NA),Beta=as.numeric(NA))[rep(NA,replications),]
	
	#Make a sequence from 0 to 1 consisting of uniformly distributed elements. These are the latent regression coefficients in our tests.
	
	betas<-unlist((sequence(replications)-1)/((replications-1)*2))
	
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
			factorloading<-runif(1)*(factorloadings[2]-factorloadings[1])+factorloadings[1]
			indicators<-data.frame(indicators,A*factorloading+rnorm(sample)*sqrt(1-factorloading^2))
			names(indicators)[[i]]<-paste("a",i,sep="")
		}
		for ( i in 1:indicatorcount){
			
			factorloading<-runif(1)*(factorloadings[2]-factorloadings[1])+factorloadings[1]
			indicators<-data.frame(indicators,B*factorloading+rnorm(sample)*sqrt(1-factorloading^2))
			names(indicators)[[indicatorcount+i]]<-paste("b",i,sep="")
	
		}

		# Standardize the indicators
		
		for (i in 1:(indicatorcount*2)){
			indicators[i]<-Make.Z(indicators[i])
		}
		
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
		
		#Calculate summed scales and run regression
	
		sumscaledata<-data.frame(A=apply(indicators[1:indicatorcount],1,mean),B=apply(indicators[(indicatorcount+1):(indicatorcount*2)],1,mean))
		
		
		regressionss<-lm.beta(lm(B~A,data=sumscaledata))

		#Calculate principal components and run regression
		pcA<-prcomp(indicators[1:indicatorcount], center=TRUE, scale.=TRUE)
		pcB<-prcomp(indicators[(indicatorcount+1):(indicatorcount*2)], center=TRUE, scale.=TRUE)
		princompdata<-data.frame(A=pcA$x[,1]*mean(pcA$rotation[,1])>0,B=pcB$x[,1]*mean(pcB$rotation[,1])>0)
		
		regressionpc<-lm.beta(lm(B~A,data=princompdata))

			
		factordata<-data.frame(A=factanal(indicators[,1:indicatorcount],1,scores="regression")$scores[,1],
		B=factanal(indicators[,(indicatorcount+1):(indicatorcount*2)],1,scores="regression")$scores[,1])
		
		regressionfa<-lm.beta(lm(B~A,data=factordata))


		#Append results to data
		results$PLS[[replication]]<-pls$path.coefs[2]
		results$SumScale[[replication]]<-regressionss[[1]]
		results$PrinComp[[replication]]<-regressionpc[[1]]
		results$Factor[[replication]]<-regressionfa[[1]]
		results$Beta[[replication]]<-betas[[replication]]
	
	}
	
	
	attach(results)
	
	pdf(file="results/figure2.pdf",width=10)
	
	par(mfrow=c(1,3)) 
	
	plot(PLS ~ Beta, asp=1, ylim=c(-.4,.6), xlim=c(0,.5))
	abline(0, 1)
	plot(SumScale ~ Beta, asp=1, ylim=c(-.4,.6), xlim=c(0,.5),ylab="Summed scales and regression")
	abline(0, 1)
	plot(Factor ~ Beta, asp=1, ylim=c(-.4,.6), xlim=c(0,.5),ylab="First common factor score and regression")
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
	
	print("start reading data")
	#library("sqldf")
	#constructData <- read.csv.sql("constructs.csv",sep="\t",row.names=FALSE)
	constructData <- read.delim("constructs.csv")
	print("read data")
	
	#Only keep the rows that exist for all of the analysis types
	
	keep<-constructData[constructData$analysis==4,c("designNumber", "replication","construct")]
	
	constructData<-merge(constructData,keep)
	
	keep<-constructData[constructData$analysis==3,c("designNumber", "replication","construct")]
	
	constructData<-merge(constructData,keep)
	rm(keep)
	
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
	

	constructData$trueScoreR2<-constructData$trueScoreCorrelation^2

	constructData$deltaR2<-constructData$deltaR2errors+constructData$deltaR2constructs

	constructData$totalR2<-constructData$trueScoreR2+constructData$deltaR2errors+constructData$deltaR2constructs
	
	#Merge the design parameters
	
	#constructData<-merge(constructData,designMatrix,by=c("designNumber"))
	
}

#
# Two plots that show the 1) cumulative function of reliability and 2) the
# density function of difference in reliabilities
#

if(!file.exists("results/figure3.pdf")){

	pdf(file="results/figure3.pdf",width=10)
	
	par(mfrow=c(1,2))

	tpls<-constructData[constructData$analysis==4,]$trueScoreR2
	tpls<-sort(tpls)[1:length(tpls) %% 1000 == 1| 1:length(tpls) == length(tpls)]
	tss<-constructData[constructData$analysis==1,]$trueScoreR2
	tss<-sort(tss)[1:length(tss) %% 100 == 1| 1:length(tss) == length(tss)]
	tf<-constructData[constructData$analysis==3,]$trueScoreR2
	tf<-sort(tf)[1:length(tf) %% 100 == 1| 1:length(tf) == length(tf)]
	
	print("1")

#	plot(1:length(tpls)/length(tpls),tpls,type="l",ylab="Share of variance explained by construct",xlab="Cumulative amount of observations")
#	lines(1:length(tss)/length(tss),tss,lty=2,type="l")
#	lines(1:length(tf)/length(tf),tf,lty=3,type="l")

	plot(density(tpls),type="l",ylab="Share of variance explained by construct",xlab="Cumulative amount of observations")
	lines(density(tss),lty=2,type="l")
	lines(density(tf),lty=3,type="l")

	minor.tick( )
	grid()
	legend(0,1, c("PLS","Summed scales","Factor scores"),lty=c(1,2,3))
	
	#Plot the difference between reliabilities of summed scale and PLS and factor scores and PLS

	tempData<-constructData[,c("designNumber","construct","replication","trueScoreR2","analysis")]

	tempData<-merge(tempData[tempData$analysis==4,], merge(tempData[tempData$analysis==1,], tempData[tempData$analysis==3,], by=c("designNumber","construct","replication")), by=c("designNumber","construct","replication"))

	print("2")

	tss<-tempData$trueScoreR2-tempData$trueScoreR2.x
	tss<-sort(tss)[1:length(tss) %% 1000 == 1 | 1:length(tss) == length(tss)]
	tf<-tempData$trueScoreR2-tempData$trueScoreR2.y
	tf<-sort(tf)[1:length(tf) %% 1000 == 1| 1:length(tf) == length(tf)]
	
	#plot(1:length(tss)/length(tss),tss,type="l",ylab="Difference in variance explained by construct",xlab="Cumulative amount of observations",lty=2)

	plot(density(tss),type="l",ylab="Difference in variance explained by construct",xlab="Cumulative amount of observations",lty=2)
	minor.tick( )
	grid()
	
	#lines(1:length(tf)/length(tf),tf,lty=3,type="l")
	lines(density(tf),lty=3,type="l")

	legend(0,.2, c("PLS vs. summed scales","PLS vs. factor scores"),lty=c(2,3))


	dev.off()
	
	# How large percent of construct scores is more affected by measurement error than the actual constructs

	rm(tempData)
}

######## TABLE 2 ############
#
# Comparing the R2s and stability of measurement
#

print("Table 2")

if(!file.exists("results/table2_full.tex")){
	
	constructDataCopy<-constructData[,c("designNumber","analysis","trueScoreR2","deltaR2constructs","deltaR2errors","totalR2","sdByData","sdByModels")]
	constructDataCopy$totalR2<-constructData$deltaR2+constructData$trueScoreR2
	
	#Eliminate seriously over fitted data
	
	predictors<-(constructData$incomingPathsCorrect+constructData$incomingPathsExtra+constructData$outgoingPathsCorrect+constructData$outgoingPathsExtra)*indicatorCounts[constructData$indicatorCount]

	include<-(predictors*5)<sampleSizes[constructData$sampleSize] & constructDataCopy$totalR2 < 1
	
	constructDataCopy[include,c("deltaR2errors","totalR2")]<-NA
	writeAlternativeComparisonTable(constructDataCopy,variables=c("trueScoreR2","deltaR2constructs","deltaR2errors","totalR2","sdByData","sdByModels"),file="table2",analysisTypes=c(4,1,3),labels=labels)
	
	rm(constructDataCopy)
}

######## TABLE 3 ############
#
#
#
#

print("Table 3")

if(!file.exists("results/table3_full.tex")){
	
	file<-"table3"
	
	
	tempData<-aggregate(constructData[,c("designNumber","analysis","trueScoreR2")], by=list(constructData$designNumber,constructData$analysis),  FUN=mean, na.rm=TRUE)
	

	tempData<-merge(tempData[tempData$analysis==4,], merge(tempData[tempData$analysis==1,], tempData[tempData$analysis==3,], by=c("designNumber")), by=c("designNumber"))

	tempData<-merge(tempData,designMatrix,by=c("designNumber"))
	
	tempData$betterThanSS<-tempData$trueScoreR2-tempData$trueScoreR2.x
	tempData$betterThanFactor<-tempData$trueScoreR2-tempData$trueScoreR2.y

	tableData<-NULL
	
	#The last column in design matrix is the design number
	
	for( i in 1:(ncol(designMatrix)-1)){
		
		varname<-colnames(designMatrix)[i]
		tableRow<-NULL
		
		for(j in 1:3){
			tableRow<-c(tableRow,mean(tempData[tempData[,varname]==j,"betterThanSS"]))
		}
		for(j in 1:3){
			tableRow<-c(tableRow,mean(tempData[tempData[,varname]==j,"betterThanFactor"]))
		}
		tableRow<-tableRow-c(rep(median(tableRow[1:3]),3),rep(median(tableRow[4:6]),3))
		tableData<-rbind(tableData,data.frame(t(tableRow),row.names=labels[[varname]]))
		
		print(tableData)
	}
	
	print(xtable(tableData,digits=3),file=paste("results/",file,"_full.tex",sep=""))
	print(xtable(tableData,digits=3),file=paste("results/",file,"_body.tex",sep=""),	hline.after=NULL,only.contents=TRUE,include.colnames=FALSE)
}


#
# Which interactions are sufficient to make PLS a better alternative than the
# alternatives. None of the first, second, third, or fourth order factorial 
# interaction determine if PLS is better.
#


if(FALSE){

	regressionData<-constructData[,c("replication", "designNumber", "construct",  "analysis","trueScoreR2")]

	regressionData<-merge(regressionData[regressionData$analysis==4,], merge(regressionData[regressionData$analysis==1,], regressionData[regressionData$analysis==3,], by=c("replication", "designNumber", "construct")), by=c("replication", "designNumber", "construct"))

	regressionData<-merge(regressionData,designMatrix,by=c("designNumber"))
	
	regressionData$betterThanSS<-regressionData$trueScoreR2-regressionData$trueScoreR2.x
	regressionData$betterThanFactor<-regressionData$trueScoreR2-regressionData$trueScoreR2.y

	regressionData$betterThanBoth<-regressionData$trueScoreR2-max(regressionData$trueScoreR2.x,regressionData$trueScoreR2.y)

	regressionData<-aggregate(regressionData,by=list(regressionData$designNumber),FUN=mean, na.rm=TRUE)
	
	testf(regressionData,"betterThanSS","betterThanSS",0,4,1)
	testf(regressionData,"betterThanFactor","betterThanFactor",0,4,1)
	testf(regressionData,"betterThanBoth","betterThanBoth",0,4,1)
	
}


######## Figure 4 ############

print("Figure 4")

if(!file.exists("results/figure4.pdf")){

#	testData<-constructData[,c("replication", "designNumber", "construct",  "analysis","trueScoreR2","meanFactorLoading", "AVE", "CR")]

#	testData<-merge(testData,modelData[,c("designNumber","replication","analysis","SRMR")],by=c("designNumber","replication","analysis"))
	
#	testData<-merge(testData[testData$analysis==4,], merge(testData[testData$analysis==1,], testData[testData$analysis==3,], by=c("replication", "designNumber", "construct")), by=c("replication", "designNumber", "construct"))

	
	# When is PLS better than summed scales or factor. 
	
	statistics<-c("meanFactorLoading","AVE","CR","SRMR")
	
	
	pdf(file="results/figure4.pdf",width=10)
	
	par(mfrow=c(2,4))
	
	#Relative or absolute
	for(i in 1:2){
		
		#Which statistic
		for(j in 1:length(statistics)){
		
			#Sum scale or factor
			for(k in 1:2){
			
				print(paste(i,j,k))
				print(statistics[j])
				testData$yvar<-testData$trueScoreR2>testData[,paste("trueScoreR2",c("x","y")[k],sep=".")]
				
				if(i==2){
					
					print("relative")
					testData$xvar<-testData[,statistics[j]]-testData[,paste(statistics[j],c("x","y")[k],sep=".")]
				}
				else{
					testData$xvar<-testData[,statistics[j]]
				}
				
				mi<-min(testData$xvar)
				ma<-max(testData$xvar)
				testData$group<-round((testData$xvar-mi)/(ma-mi)*1000)
				
				print(summary(testData[c("xvar","yvar")]))
				graphData<-aggregate(testData[,c("xvar","yvar")],by=list(testData$group),FUN=mean, na.rm=TRUE)

				print(summary(graphData[c("xvar","yvar")]))
				
				smoothed<-loess.smooth(graphData$xvar,graphData$yvar, family = "gaussian",span=.25)
				
				if(k==1){
					templabel<-statistics[j]
					if(j==1){
						templabel<-"Mean factor loading"
					}
					
					if(i==1){
						templabel<-paste(templabel,", absolute",sep="")
					}
					else{
						templabel<-paste(templabel,", difference",sep="")

					}
					plot(smoothed$x,smoothed$y,type="l",xlab=templabel,ylab="Likelihood that PLS score is more reliable",ylim=c(0,1),lty=2)
					grid()
				}
				else{
					lines(smoothed$x,smoothed$y,lty=3,type="l")
					if(j==1 & i==1){
						legend(0,1, c("PLS vs. summed scales","PLS vs. factor scores"),lty=c(2,3))
					}
				}
			}
		}
	}
	dev.off()
}

if(!file.exists("results/table5_full.tex") & FALSE){

#	testData<-constructData[,c("replication", "designNumber", "construct",  "analysis","trueScoreR2","meanFactorLoading", "AVE", "CR")]

#	testData<-merge(testData,modelData[,c("designNumber","replication","analysis","SRMR")],by=c("designNumber","replication","analysis"))
	
#	testData<-merge(testData[testData$analysis==4,], merge(testData[testData$analysis==1,], testData[testData$analysis==3,], by=c("replication", "designNumber", "construct")), by=c("replication", "designNumber", "construct"))

	
	# When is PLS better than summed scales or factor. 
	
	statistics<-c("meanFactorLoading","AVE","CR","SRMR")
	
	
	#Which statistic
	for(j in 1:length(statistics)){
	
		#Sum scale or factor
		for(k in 1:2){
		
			print(paste(i,j,k))
			print(statistics[j])
			testData$yvar<-testData$trueScoreR2>testData[,paste("trueScoreR2",c("x","y")[k],sep=".")]
			
			if(i==2){
				
				print("relative")
				testData$xvar<-testData[,statistics[j]]-testData[,paste(statistics[j],c("x","y")[k],sep=".")]
			}
			else{
				testData$xvar<-testData[,statistics[j]]
			}
			
			mi<-min(testData$xvar)
			ma<-max(testData$xvar)
			testData$group<-round((testData$xvar-mi)/(ma-mi)*1000)
			
			print(summary(testData[c("xvar","yvar")]))
			graphData<-aggregate(testData[,c("xvar","yvar")],by=list(testData$group),FUN=mean, na.rm=TRUE)

			print(summary(graphData[c("xvar","yvar")]))
			
			smoothed<-loess.smooth(graphData$xvar,graphData$yvar, family = "gaussian",span=.25)
			
			if(k==1){
				templabel<-statistics[j]
				if(j==1){
					templabel<-"Mean factor loading"
				}
				
				if(i==1){
					templabel<-paste(templabel,", absolute",sep="")
				}
				else{
					templabel<-paste(templabel,", difference",sep="")

				}
				plot(smoothed$x,smoothed$y,type="l",xlab=templabel,ylab="Likelihood that PLS score is more reliable",ylim=c(0,1),lty=2)
				grid()
			}
			else{
				lines(smoothed$x,smoothed$y,lty=3,type="l")
				if(j==1 & i==1){
					legend(0,1, c("PLS vs. summed scales","PLS vs. factor scores"),lty=c(2,3))
				}
			}
		}
		
	}
	dev.off()
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
	
	# This is not strictly speaking the right way to determine if a path is 
	# speficied because also non-convergent analyses result in NA estimates. 
	# A more approriate way would be to read the tested model from MapReduce 
	# input data or to test if the estimate is NA for all analysis types.
	

	relationshipData$specifiedAsPath<-!is.na(relationshipData$regressionEstimate)
	

	
}

if( ! exists("oldRelationshipData")){
	oldRelationshipData <- read.delim("olddata2/relationships.csv")
	oldRelationshipData$regressionARE<-abs(oldRelationshipData$regressionTrueScore-oldRelationshipData$regressionEstimate)

	oldRelationshipData$correlationBias<-oldRelationshipData$estimatedCorrelation-oldRelationshipData$trueCorrelation*oldRelationshipData$correlationAttenuationCoefficient
	oldRelationshipData$correlationError<-abs(oldRelationshipData$estimatedCorrelation-oldRelationshipData$trueCorrelation)

	tempMatrix<-data.frame(designNumber=1:729,designMatrix)
	oldRelationshipData<-merge(oldRelationshipData,tempMatrix[,c("designNumber","sampleSize")])
	
	oldRelationshipData$sampleSize<-sampleSizes[oldRelationshipData$sampleSize]
	oldRelationshipData$t<-oldRelationshipData$regressionEstimate/oldRelationshipData$regressionSE
	oldRelationshipData$p<-(1-pt(abs(oldRelationshipData$t),oldRelationshipData$sampleSize-1))*2
	
	# This is not strictly speaking the right way to determine if a path is 
	# speficied because also non-convergent analyses result in NA estimates. 
	# A more approriate way would be to read the tested model from MapReduce 
	# input data or to test if the estimate is NA for all analysis types.
	

	oldRelationshipData$specifiedAsPath<-!is.na(oldRelationshipData$regressionEstimate)
	

	
}

#
# Two plots that show the 1) cumulative function of reliability and 2) the
# density function of difference in reliabilities
#

print("Figure 5")

if(!file.exists("results/figure5.pdf")){

	pdf(file="results/figure5.pdf",width=10)
	
	par(mfrow=c(1,2))

	tpls<-relationshipData[relationshipData$analysis==4 & ! is.na(relationshipData$regressionARE),]$regressionARE
	tpls<-sort(tpls)[1:length(tpls) %% 1000 == 1  | 1:length(tpls) ==length(tpls) ]
	tss<-relationshipData[relationshipData$analysis==1 & ! is.na(relationshipData$regressionARE),]$regressionARE
	tss<-sort(tss)[1:length(tss) %% 1000 == 1  | 1:length(tss) ==length(tss) ]
	tf<-relationshipData[relationshipData$analysis==3 & ! is.na(relationshipData$regressionARE),]$regressionARE
	tf<-sort(tf)[1:length(tf) %% 1000 == 1 | 1:length(tf) ==length(tf) ]
	
	print("1")

	plot(1:length(tpls)/length(tpls),tpls,type="l",ylab="Absolute relative error",xlab="Cumulative amount of observations",ylim=c(0,1))
	minor.tick( )
	grid()
	lines(1:length(tss)/length(tss),tss,lty=2,type="l")
	lines(1:length(tf)/length(tf),tf,lty=3,type="l")
	legend(0,1, c("PLS","Summed scales","Factor scores"),lty=c(1,2,3))
	
	#Plot the difference between reliabilities of summed scale and PLS and factor scores and PLS

	tempData<-relationshipData[,c("designNumber","to","from","replication","regressionARE","analysis")]

	tempData<-merge(tempData[tempData$analysis==4,], merge(tempData[tempData$analysis==1,], tempData[tempData$analysis==3,], by=c("designNumber","to","from","replication")), by=c("designNumber","to","from","replication"))

	print("2")

	tss<-tempData$regressionARE-tempData$regressionARE.x
	tss<-tss[!is.na(tss)]
	tss<-sort(tss)[1:length(tss) %% 1000 == 1]
	tf<-tempData$regressionARE-tempData$regressionARE.y
	tf<-tf[!is.na(tf)]
	tf<-sort(tf)[1:length(tf) %% 1000 == 1]
	
	plot(1:length(tss)/length(tss),tss,type="l",ylab="Difference in absolute relative error",xlab="Cumulative amount of observations",lty=2)
	minor.tick( )
	grid()
	
	lines(1:length(tf)/length(tf),tf,lty=3,type="l")
	legend(0,.75, c("PLS vs. summed scales","PLS vs. factor scores"),lty=c(2,3))


	dev.off()
	
	rm(tempData)
}

######## TABLE 2 ############
#
# Comparing the R2s and stability of measurement
#

print("Table 4")

if(!file.exists("results/table4_full.tex")){
	
	
	# Type I error is when the population regression coefficient is either zero 
	# or of different sign than the estimated coefficient
	
	oldRelationshipData$TypeI05<-oldRelationshipData$p<=.05 & (oldRelationshipData$regressionEstimate*oldRelationshipData$regressionTrueScore <=0)
	
	oldRelationshipData$TypeI01<-oldRelationshipData$p<=.01 & (oldRelationshipData$regressionEstimate*oldRelationshipData$regressionTrueScore <=0)

	oldRelationshipData$TypeI001<-oldRelationshipData$p<=.001 & (oldRelationshipData$regressionEstimate*oldRelationshipData$regressionTrueScore <=0)

	# Type II error is when the population regression coefficient is different from zero but the p-value is non-significant sign than the estimated coefficient

	oldRelationshipData$TypeII05<-oldRelationshipData$p>.05 & oldRelationshipData$regressionTrueScore != 0

	oldRelationshipData$TypeII01<-oldRelationshipData$p>.01 & oldRelationshipData$regressionTrueScore != 0

	oldRelationshipData$TypeII001<-oldRelationshipData$p>.001 & oldRelationshipData$regressionTrueScore != 0

	print("Prepared data")

	writeAlternativeComparisonTable(oldRelationshipData,variables=c("regressionARE","TypeI01","TypeII01"),file="table4",analysisTypes=c(5,1,3),labels=labels)
	
	
}
