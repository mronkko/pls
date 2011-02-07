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

labels<-list(CR="CR",AVE="Root AVE",minFactorLoading="Minumum factor loading",meanFactorLoading="Mean factor loading",
	maxCrossLoading="Maximum cross loading",AVEMinusMaxCorrelation="Root AVE - max correlation",sumscale="Summed scales",component="Components",factor="Factors",pls="PLS")

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
	
	tempData$trueScoreCorrelation<-abs(tempData$trueScoreCorrelation)
	tempData$minFactorLoading<-abs(tempData$minFactorLoading)
	tempData$meanFactorLoading<-abs(tempData$meanFactorLoading)
	tempData$maxCrossLoading<-abs(tempData$maxCrossLoading)

	# We typically examine the square root of AVE

	tempData$AVE<-sqrt(tempData$AVE)
	tempData$AVEMinusMaxCorrelation<-tempData$AVE-tempData$maxCorrelationWithOtherConstruct

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

if(!file.exists("results/table8_full.tex")){

	dependents<-c("deltaR2","sdByData","correlationBias","correlationError","regressionARE","regressionSE")
	
	allVars<-c("analysis","designNumber",dependents)
	
	# Rare events logistic regression on when PLS is better than others.
	
	tempdata<-aggregate(constructData[,intersect(names(constructData),allVars)], by=list(constructData$designNumber,constructData$analysis),  FUN=mean, na.rm=TRUE)
	
	tempdata<-merge(tempdata,aggregate(relationshipData[,intersect(names(relationshipData),allVars)], by=list(relationshipData$designNumber,relationshipData$analysis),  FUN=mean, na.rm=TRUE))

	tempdata<-reshape(tempdata,v.names=variables,idvar="designNumber",timevar="analysis",direction="wide")

	#Order the data bu design number

	tempdata<-tempdata[order(tempdata[,"designNumber"]),]

	
	tempdata2<-data.frame(tempdata[,"designNumber"])
	
	for(i in 1:length(dependents)){
		tempdata2[,dependents[i]]<-tempdata[,paste(dependents[i],4,sep=".")]==apply(tempdata[,paste(dependents[i],1:4,sep=".")],1,min)
	}

	#Merge the design numbers
	
	tempdata2<-cbind(tempdata2,designMatrix)

	relogit<-list()
	for(i in 1:length(dependents)){
		
		relogit[[length(relogit)+1]]<-zelig(as.formula(paste(dependents[i],"~",paste(colnames(designMatrix),collapse=" + "))), model="relogit",data=tempdata2,tau=c(mean(tempdata2[,dependents[i]])))
		
		# TODO: relogit postestimation
		
	}
	#Very cumbersome, but works.
	for(i in 1:6) assign(paste("m",i,sep=""),relogit[[i]])
	write(format(mtable(m1,m2,m3,m4,m5,m6),forLaTeX=TRUE),file="results/table8_full.tex")
}

######## TABLE 9 ############

if(FALSE & !file.exists("results/table9_full.tex")){

	# Regression on all quality measures. Not grouped by experimental 
	# conditions. This results in a huge regression table
	
	# Precision and accuracy
	tempData<-relationshipData[,c("replication","designNumber","to","from","analysis","correlationAttenuationCoefficient","regressionTrueScore","regressionEstimate","regressionSE")]
	tempData$regressionARE<-abs(tempData$regressionTrueScore-tempData$regressionEstimate)
	writeComparisonTable(tempData,variables=c("regressionARE","regressionSE"),file="table9",analysisTypes=analysisTypes,labels=labels)
}

####### DISTRIBUTION PLOTS ########

#
# Show that the distribution of parameter estimates around the true value do not # follow the t-distribution. Use only data without method variance here and
# include correctly specified models. Only data with 100 observations
#


if(FALSE){
	CorretModelsWithoutMethodVariance<-relationshipData[relationshipData$designNumber %in% which((designMatrix[,4]==1) & (designMatrix[,5]==1) & (designMatrix[,11]==1)&(designMatrix[,6]==2)),]
	
	# Do a four by four matrix of distribution plots
	
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



















