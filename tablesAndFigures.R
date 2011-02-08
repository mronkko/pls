#
# This file creates a set of simulated data and runs PLS, SEM, and SummedScales on these.
#


library(plspm)
library(lattice)
library(lavaan)
library(car)
library(Zelig)

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
	pls="PLS",
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
	methodVariance="Method variances")

designMatrix<-createdesignMatrix()


#Only draw the tables and figures that do not yet exist.

# Figures 1, 2, and 3 are drawn with PowerPoint

######## FIGURE 4 ############

if(!file.exists("results/figure4.pdf")){
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

######## TABLE 3 ############
if(!file.exists("results/table3_full.tex")){
	# General information about quality of measurement
	
	writeDescriptivesTable(constructData,variables=c("CR","AVE","minFactorLoading","meanFactorLoading",
	"maxCrossLoading","AVEMinusMaxCorrelation"),file="table3",analysisTypes=analysisTypes,labels=labels)
}

######## TABLE 4 ############

if(!file.exists("results/table4_full.tex")){

	# Construct score reliablity and validity (Hypothesis 1)
	

	writeComparisonTable(constructData,variables=c("trueScoreCorrelation","deltaR2"),file="table4",analysisTypes=analysisTypes,labels=labels)
	
	# Show the experimental conditions in which PLS was best
}

######## TABLE 5 ############

if(!file.exists("results/table5_full.tex")){

	# Construct score stability (Hypotheses 1)
	
	writeComparisonTable(constructData,variables=c("sdByData","sdByModels"),file="table5",analysisTypes=analysisTypes,labels=labels)
}

#Only read in this data if it is not in memory already. It takes a while to read

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
}


######## TABLE 6 ############

if(!file.exists("results/table6_full.tex")){

	# Correlations (Hypothesis 2)
	
	# Attenuation and correlationBias

	writeComparisonTable(relationshipData,variables=c("correlationAttenuationCoefficient","correlationBias","correlationError"),file="table6",analysisTypes=analysisTypes,labels=labels)
}

######## TABLE 7 ############

if(!file.exists("results/table7_full.tex")){

	# Regression coefficients (Hypothesis 3)
	
	# Precision and accuracy
	writeComparisonTable(relationshipData,variables=c("regressionARE","regressionSE"),file="table7",analysisTypes=analysisTypes,labels=labels)
}

######## TABLE 8 ############

if( !file.exists("results/table8_full.tex")){

	dependents<-c("deltaR2","sdByData","correlationBias","correlationError","regressionARE","regressionSE")
	
	allVars<-c("analysis","designNumber",dependents)
	
	# Crosstabulate 
	
	tempdata<-aggregate(constructData[,intersect(names(constructData),allVars)], by=list(constructData$designNumber,constructData$analysis),  FUN=mean, na.rm=TRUE)
	
	tempdata<-merge(tempdata,aggregate(relationshipData[,intersect(names(relationshipData),allVars)], by=list(relationshipData$designNumber,relationshipData$analysis),  FUN=mean, na.rm=TRUE))

	tempdata<-reshape(tempdata,v.names=variables,idvar="designNumber",timevar="analysis",direction="wide")

	#Order the data bu design number

	tempdata<-tempdata[order(tempdata[,"designNumber"]),]

	
	tempdata2<-data.frame(tempdata[,"designNumber"])
	
	for(i in 1:length(dependents)){
		tempdata2[,dependents[i]]<-tempdata[,paste(dependents[i],4,sep=".")]==apply(tempdata[,paste(dependents[i],1:4,sep=".")],1,min)
	}

	# Merge the design numbers
	

	# Loop over dependents and generate vectors of frequencies
	
	resultTable<-NULL
	for(i in 1:length(dependents)){
		resultCol<-NULL
		temp<-designMatrix[tempdata2[,dependents[i]]==TRUE,]
		for(i in 1:ncol(temp)){
			for(k in 1:3){
				resultCol<-c(resultCol,sum(temp[,i]==k)/243)
			}
		}
		resultTable<-cbind(resultTable,matrix(resultCol,ncol=1))
		
	}
	
	#Add row names 
	rownames=c(numberOfConstructs,expectedNumberOfOutgoingPaths,paste(populationPathValues),paste(omittedPathsShare*100,"%",sep=""),extraPaths,sampleSizes,indicatorCounts,factorLoadings,factorLoadingIntervals,maxErrorCorrelations,methodVariances)
	
	tableData<-data.frame(rownames,resultTable)
	print(resultTable)

	command<-NULL
	pos<-NULL
	for(i in 1:ncol(designMatrix)){
		pos=c(pos,(i-1)*3)
		command=c(command,paste("\\multicolumn{",length(dependents)+1,"}{l}{",labels[[colnames(designMatrix)[i]]],"}\\\\"))

	}
	
	file="table8"
	add.to.row=list(as.list(pos),command)
	print(xtable(tableData,digits=3),file=paste("results/",file,"_full.tex",sep=""),add.to.row=add.to.row)
	print(xtable(tableData,digits=3),file=paste("results/",file,"_body.tex",sep=""),include.rownames=FALSE,
hline.after=NULL,only.contents=TRUE,include.colnames=FALSE,add.to.row=add.to.row)
}

######## TABLES 9 ############

if(FALSE & !file.exists("results/table9_full.tex")){

	# Correlation table where conditions are on columns and all construct related things are on rows
	
}

######## TABLES 10 ############

if(FALSE & !file.exists("results/table10_full.tex")){
	# Two-level regression table where all construct related things are 
	# dependents and all model related things are on the second level and 
	# all construct related things on the first level

}

######## TABLES 11 ############

if(FALSE & !file.exists("results/table11_full.tex")){
	# Correlation table where all relationship related things are included
}

######## TABLES 12 ############

if(FALSE & !file.exists("results/table12_full.tex")){
	# Two-level regression table where all construct related things are 
	# dependents and all model related things are on the second level and 
	# all construct related things on the first level

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
	for(i in 1:length(analysisTypes)){

		tempData <- dataSample[dataSample$analysis==i,]

		plot(tempData$trueR2,tempData$estimatedR2,xlab="Real R2",ylab="Estimated R2",main=labels[[analysisTypes[i]]])
		abline(0, 1)

	}
	dev.off()
}

if(!file.exists("results/figure6.pdf")){
	
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
if(!file.exists("results/figure7.pdf")){
	
	# Do a distribution plot for p-value when there is no effect
	# Also there must be only one incoming or outgoing path.
	
	#pdf(file="results/figure6.pdf")
	
	#Choose the construct where there is only one path to other construct
	tempConstructs<-cbind(constructData[constructData$incomingPathsExtra+constructData$incomingPathsCorrect+constructData$outgoingPathsExtra+constructData$outgoingPathsCorrect==1,c("construct","replication","designNumber","incomingPathsExtra","incomingPathsCorrect","outgoingPathsExtra","outgoingPathsCorrect")],TRUE)
	
	
	#Only use data where the relationship is very close to zero
	dataSample<-relationshipData[abs(relationshipData$trueCorrelation)<.01,]
	
	tempMatrix<-data.frame(designNumber=1:729,designMatrix)
	
	useDesigns<-tempMatrix[tempMatrix$methodVariance==1,]$designNumber

	# Only use data without method variance
	
	dataSample<-dataSample[dataSample$designNumber %in% useDesigns,]
	
	# Remove NAs
	dataSample<-dataSample[!is.na(dataSample$p),]

	dataSample<-merge(dataSample,tempConstructs,by.y=c("construct","replication","designNumber"),by.x=c("to","replication","designNumber"),all.x=TRUE)
	dataSample<-merge(dataSample,tempConstructs,by.y=c("construct","replication","designNumber"),by.x=c("from","replication","designNumber"),all.x=TRUE)
	
	dataSample<-dataSample[dataSample$TRUE.x==TRUE | dataSample$TRUE.y==TRUE,]
	
	print(names(dataSample))
	print(dataSample[1:100,c("trueCorrelation","regressionEstimate","regressionSE","t","p","incomingPathsExtra.x","incomingPathsCorrect.x","outgoingPathsExtra.x","outgoingPathsCorrect.x","incomingPathsExtra.y","incomingPathsCorrect.y","outgoingPathsExtra.y","outgoingPathsCorrect.y")],digits=3)

	stop("debug")

	par(mfrow=c(2,2)) 
	
	for(i in 1:length(analysisTypes)){

		tempData <- as.vector(dataSample[dataSample$analysis==i,"p"])
		tempData <- tempData[order(tempData)]
		
		x<-((1:length(tempData)))/(length(tempData))

		plot(x,tempData,xlab="Cumulative propability",ylab="Estimated value",main=labels[[analysisTypes[i]]],log="xy")
		abline(0,1)

	}
	#dev.off()
}
if(FALSE){
	CorretModelsWithoutMethodVariance<-relationshipData[relationshipData$designNumber %in% which((designMatrix[,4]==1) & (designMatrix[,5]==1) & (designMatrix[,11]==1)&(designMatrix[,6]==2)),]
	
	# Do four of distribution plots
	
	#pdf(file="results/figure1.pdf")
	
	limits<-c(.1,.25,.4)
	
	par(mfrow=c(length(analysisTypes),length(limits)+1)) 
	
	for(j in 1:(length(limits)+1)){
		# Rows
		dataSample<-CorretModelsWithoutMethodVariance
		if(j<=length(limits)){
			dataSample<-dataSample[abs(dataSample$regressionTrueScore)<=limits[j],]
		}
		if(j>1){
			dataSample<-dataSample[abs(dataSample$regressionTrueScore)>limits[j-1],]
		}
		# Columns
		for(i in 1:length(analysisTypes)){
			d <- density((dataSample$regressionTrueScore-dataSample$regressionEstimate)[dataSample$analysis==i],na.rm=TRUE)
			plot(d)
		}
	}
	#dev.off()
}

#
# PRAGMATIC PART: TYPE I AND TYPE II ERROR RATE
#
#



















