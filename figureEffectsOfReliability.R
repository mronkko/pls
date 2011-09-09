#! /usr/bin/env Rscript
#
# This file creates a set of simulated data and runs PLS and SummedScales on
# these. The results are then used to draw a figure that illustrates the effect # of reliability on the parameter estimates
#


#library(car)
library(QuantPsyc)
library(Hmisc)
library(multicore)

# Read simulation parameters
source("include/parameters.R")

# Read function definitions
source("include/functionsPrepare.R")
source("include/functions.R")


# Start by generating the simulation samples

#
# We start by generating a fractional factorial design matrix. Running this 
# simulation with full factorial design is not feasible due to the number of 
# factors.
#


designMatrix<-createdesignMatrix()

replications<-100
maxThreads<-20

doOneReplication<-function(){

    # Generate all population models and all tested models in mapper and send these to reducer that will generate the construct and indicator scores and test the models.

	# There are 27 different parametrizations to the population model. The
	# Parametrization changes every 27 models
	
	
	#Object for storing the results
	dataobj<-NULL
	
	for(design in 1:729){
	
		# Use only the first values of factors thar are not modeled in this simulation.

		thisDesignRow<-designMatrix[design,]
		
    	if(thisDesignRow["populationPathValues"]!=1 || thisDesignRow["maxErrorCorrelation"]!=1 || thisDesignRow["methodVariance"]!=1)  next 
    	
		
		constructCount=numberOfConstructs[thisDesignRow["numberOfConstructs"]]
		indicatorCount=indicatorCounts[thisDesignRow["indicatorCount"]]
		sampleSize=sampleSizes[thisDesignRow["sampleSize"]]

		# Data generation fails if we try to make too few observations
		
		if(sampleSize<=constructCount*indicatorCount) next
		
    	#print(paste("Not skipping design",design))
    	

    	# Set the population parameters
    	thisNumberOfConstructs<-numberOfConstructs[thisDesignRow[1]]
    	thisExpectedNumberOfOutgoingPaths<-expectedNumberOfOutgoingPaths[thisDesignRow[2]]
    	

		# Generate a population model here. We need two matrices. 
		# First matrix contains 1s and 0s to store which paths are specified
		# Second matrix contains the path values

		populationModel<-NULL
	  	
	  	while(is.null(populationModel)){
	  	
		  	populationModelWhichPaths <- generateRandomModel(thisNumberOfConstructs,thisExpectedNumberOfOutgoingPaths)

			populationModel<-setPopulationModelPathValues(populationModelWhichPaths,c(-.3,-.2,-.1,0,0,.1,.2,.3),discrete=TRUE)
		}
	    
		# Generate the tested model
		testedModel<-generateTestedModel(populationModelWhichPaths,omittedPathsShare[thisDesignRow["omittedPathsShare"]],extraPaths[thisDesignRow["extraPaths"]])

		
		
		constructTrueScores <- mvrnorm(n=sampleSizes[thisDesignRow["sampleSize"]],rep(0,constructCount),populationModel$covariances)

		seed<-.Random.seed

		indicators<-generateData(constructTrueScores,indicatorCounts[thisDesignRow[7]],factorLoadings[thisDesignRow[8]],factorLoadingIntervals[thisDesignRow[9]],0,0,uncorrelatedRandomErrors=FALSE)
		
		# Restore random number generator to previous state
		.Random.seed<-seed
		indicatorsNoErrorCorrelations<-generateData(constructTrueScores,indicatorCounts[thisDesignRow[7]],factorLoadings[thisDesignRow[8]],factorLoadingIntervals[thisDesignRow[9]],0,0,uncorrelatedRandomErrors=TRUE)

	
		#Start by constructing with inner and outer model matrixes. 
			
		outer <- list()
	
		for ( i in 1:constructCount ){
			outer[[i]] <- c(((i-1)*indicatorCount+1):(i*indicatorCount))
		}
	
		# plspm does not like NA:s in the matrix, so we will replace these with 
		# zeros.
		testedModel[is.na(testedModel)]<-0

		colnames(testedModel)<-paste("C",c(1:constructCount),sep="")
		rownames(testedModel)<-paste("C",c(1:constructCount),sep="")

		plsResults<-NULL

		tryCatch(
		
			plsResults<-plspm(indicators$indicators,testedModel,outer, rep("A",constructCount), scheme= "path")
				,error = function(e){
					debugPrint(e)
				}
			)


		if(! is.null(plsResults)){
		
			tempres<-list()

			# summed scales with data with uncorrelated errors
			tempres[[1]]<-estimateWithRegression(testedModel,indicatorsNoErrorCorrelations$indicators,"sumscale")

			# calculate summed scale estimates
			tempres[[2]]<-estimateWithRegression(testedModel,indicators$indicators,"sumscale")

			# summed scales with data with uncorrelated errors weighted with the PLS weights
			tempres[[3]]<-estimateWithRegression(testedModel,t(t(indicatorsNoErrorCorrelations$indicators)*plsResults$out.weights),"sumscale")
			
			reliabilities<-cor(plsResults$latents,constructTrueScores)
			
			# Calculate disattennuated correlation matrix
			
			# Start by calculating a vector of cronbach alphas
			
			alphas<-NULL
			sumscales<-NULL
			
			for(i in 1:constructCount){
				indicatorBlock<-indicators$indicators[,((i-1)*indicatorCount+2):(i*indicatorCount)]
				
				alphas<-c(alphas,alpha(indicatorBlock)$std.alpha
				sumscales<-cbind(sumscales,rowSums(indicatorBlock))
			}
			
			# Calculate correlations of summed scales
			sumscaleCorralations<-cor(sumScales)
			
			# Disattennuated correlations
			
			disattennuatedCorrelations<-sumscaleCorrelations/sqrt(alphas %o% alphas)
			
			# Calculate the tested model coefficients with disattennuated correlation coefficients
			
			for(c in 1:ncol(plsResults$path.coefs)){
				for(r in 1:nrow(plsResults$path.coefs)){
					if(plsResults$path.coefs[r,c]!=0){

						
						#Store the following data for each estimated path
						# 1 True reliability of the "from" construct
						# 2 True reliablity of the "to" construct
						# 3 True path coefficient 
						# 4 True covariance between constructs 
						# 5 Path coefficient with uncorrelated errors with regression
						# 6 Final path coefficient estimated with regression
						# 7 Path coefficient with uncorrelated errors with PLS
						# 8 Final path coefficient estimated with PLS


						datarow<-c(reliabilities[r,r],reliabilities[c,c],populationModel$paths[r,c],populationModel$covariances[r,c])

						for(i in 1:3){
							to<-tempres[[i]]$paths[,"To"]==paste("C",r,sep="")
							from<-tempres[[i]]$paths[,"From"]==paste("C",c,sep="")
							
							datarow<-c(datarow, tempres[[i]]$paths[ to & from ,"Estimate"])

						}
						datarow<-c(datarow,plsResults$path.coefs[r,c])
						
						# Calculate the regression coefficient from disattenuated correlations
						
						
						dataobj<-rbind(dataobj,datarow)
					}
				}
			}
		}
		
	}
	#print("RETURNING")
	#print(data)
	return(dataobj)
}



rm("dataobj")

if(!exists("dataobj")){


	dataobj<-NULL

	jobs<-list()
	if(FALSE){
	for(replicationNumber in 1:replications){
		print(paste("Starting replication", replicationNumber))
		
		jobs[[(replicationNumber-1)%%maxThreads+1]]<-parallel(doOneReplication(),mc.set.seed=TRUE)
		
		if(replicationNumber%%maxThreads==0){
			print("Collecting results")
			results<-collect(jobs)
			for(i in 1:length(results)){
				dataobj<-rbind(dataobj,results[[i]])
				
			}
			save(dataobj,file=paste("dataobj",replicationNumber,sep=""))
			dataobj<-NULL

			jobs<-list()
		}
	}
	}
	temp<-NULL
	
	for(replicationNumber in 1:replications){
		if(replicationNumber%%maxThreads==0){

			print("Loading results")
			load(file=paste("dataobj",replicationNumber,sep=""))

			temp<-rbind(temp,dataobj)
				
		}
	}
	dataobj<-data.frame(temp)
}

#Print graphs using the data
dataobj$reliability<-(dataobj[,1]+dataobj[,2])/2

pdf(file="results/figureX.pdf",width=10)
	
par(mfrow=c(4,4))
par(mar = c(2, 2, .5, .5)) 
for(r in 1:4){
	if(r==1) tempdata<-dataobj[dataobj$reliability<.5,]
	else if(r==2) tempdata<-dataobj[dataobj$reliability>=.5&dataobj$reliability<.7,]
	else if(r==3) tempdata<-dataobj[dataobj$reliability>=.7&dataobj$reliability<.85,]
	else if(r==4) tempdata<-dataobj[dataobj$reliability>=.85,]

	for(e in 1:4){
		thisdata<-tempdata[(abs(tempdata[,3])==(e-1)/10 ) & (abs(tempdata[,4])==(e-1)/10 ),]
		flip<-(thisdata[,3]==abs(thisdata[,3]))*2-1
		
		plot(density(thisdata[,7]*flip),lty=2,main=NA,xlab=NA,ylab=NA)
		lines(density(thisdata[,8]*flip),lty=1,type="l")
		lines(density(thisdata[,5]*flip),lty=2,type="l",col=rgb(0.7,0.7,0.7))
		lines(density(thisdata[,6]*flip),lty=1,type="l",col=rgb(0.7,0.7,0.7))
		abline(v=(e-1)/10)
		if(e==1 & r==1){
			legend("topleft",c("PLS no errors","PLS w errors","SS no errors","SS w errors"),lty=c(2,1,2,1),col=c(rgb(0,0,0),rgb(0,0,0),rgb(0.7,0.7,0.7),rgb(0.7,0.7,0.7)))
		}
	}
}
dev.off()

