#
# This file creates a set of simulated data and runs PLS, SEM, and SummedScales on these.
#


library(plspm)
library(lattice)
library(lavaan)
library(reshape)
library(car)

# Read simulation parameters
source("include/parameters.R")

# Read function definitions
source("include/functions.R")
source("include/functionsTablesAndFigures.R")

#Only draw the tables and figures that do not yet exist.


######## FIGURE 1 ############

if(!file.exists("results/figure1.pdf")){
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
	results<-data.frame(PLS=as.numeric(NA),SEM=as.numeric(NA),Regression=as.numeric(NA),RealBeta=as.numeric(NA))[rep(NA,replications),]
	
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
				,error = function(e){
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
		results$RealBeta[[replication]]<-betas[[replication]]
	
	}
	
	
	attach(results)
	
	pdf(file="results/figure1.pdf")
	
	par(mfrow=c(1,3)) 
	
	plot(PLS ~ RealBeta, asp=1, ylim=c(-.4,1.1), xlim=c(0,1))
	abline(0, 1)
	plot(SumScale ~ RealBeta, asp=1, ylim=c(-.4,1.1), xlim=c(0,1))
	abline(0, 1)
	plot(SEM ~ RealBeta, asp=1, ylim=c(-.4,1.1), xlim=c(0,1))
	abline(0, 1)
	dev.off()
}

#Only read in this data if it is not in memory already. It takes a while to read

if( ! exists("constructData")){
	constructData <- read.delim("data/constructs.csv")
}

######## TABLE 1 ############

#
# How frequently the correlation with true score is negative
#
# This does not go into the table to save space. Just print to output 

for(i in 1:length(analysisTypes)){
	count<-sum(constructData$analysis==i)
	negatives<-sum(constructData$analysis==i & constructData$trueScoreCorrelation<0) 
	print(paste("Analysis:",analysisTypes[i]," share of negative truescore correlations ",negatives,"/",count,"=",negatives/count))
}

# Fix the data so that negative correlations are not an issue
tempData<-constructData[,c("replication","designNumber","construct","analysis","trueScoreCorrelation","deltaR2")]

tempData$trueScoreCorrelation<-abs(tempData$trueScoreCorrelation)

writeComparisonTable(tempData,variables=c("trueScoreCorrelation","deltaR2"),file="table1",analysisTypes=analysisTypes)
























