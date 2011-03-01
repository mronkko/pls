#! /usr/bin/env Rscript

#
# Prepares an input file for the MapReduce jobs. The script runs several
# replications and for each replication, it generates
# one population model for each combination of population model parametrization
# (number of constructs, paths set to zero, population path values). The it
# generates all combinations of tested model parameters (paths omitted, extra
# paths) and passes the population model and tested model combinations to
# the reducer that will create datasetss based on the population model and
# evaluate the tested model with these data.
#
#

# Ensure that there is no garbage to mess things up
rm(list = ls(all = TRUE))

set.seed(Sys.time())
 
# Read simulation parameters
source("include/parameters.R")

# Read function definitions
source("include/functionsPrepare.R")
source("include/functions.R")


#
# We start by generating a fractional factorial design matrix. Running this 
# simulation with full factorial design is not feasible due to the number of 
# factors.
#


designMatrix<-createdesignMatrix()

# Keep a counters of reduce jobs
counter<-0

for(replicationNumber in 1:replications){

    # Generate all population models and all tested models in mapper and send these to reducer that will generate the construct and indicator scores and test the models.

	# There are 27 different parametrizations to the population model. The
	# Parametrization changes every 27 models
	
	for(populationModelIndex in 1:27){
    	
    	# Set the population parameters
    	thisNumberOfConstructs<-numberOfConstructs[designMatrix[populationModelIndex*27,1]]
    	thisExpectedNumberOfOutgoingPaths<-expectedNumberOfOutgoingPaths[designMatrix[populationModelIndex*27,2]]
    	
    	# populationPathValues is a list. All other parameters are stored in 
    	# vectors
    	thisPopulationPathValues<-populationPathValues[[designMatrix[populationModelIndex*27,3]]]

		# Generate a population model here. We need two matrices. 
		# First matrix contains 1s and 0s to store which paths are specified
		# Second matrix contains the path values

#		print("Choosing paths")
	  
		populationModel<-NULL
	  	
	  	while(is.null(populationModel)){
	  	
		  	populationModelWhichPaths <- generateRandomModel(thisNumberOfConstructs,thisExpectedNumberOfOutgoingPaths)

			populationModel<-setPopulationModelPathValues(populationModelWhichPaths,thisPopulationPathValues)
		}
	    
		# Generate all tested models for this population model.

		testedModels<-list()

		for(omittedPathsShareIndex in 1:3){
		    for(extraPathsIndex in 1:3){

				testedModels[[(omittedPathsShareIndex-1)*3+extraPathsIndex]] <- generateTestedModel(populationModelWhichPaths,omittedPathsShare[omittedPathsShareIndex],extraPaths[extraPathsIndex])
			}
		}
		
		debugPrint("Population model")
		debugPrint(populationModelWhichPaths)
		debugPrint("Tested model")
		debugPrint(testedModels)
		
		# The population model is tested with 3 different sample sizes. Run each
		# as a separeate reduce job.
		
		for(sampleSizeIndex in 1:3){

			counter <- counter +1


			#write a progress indicator on screen
			print(paste(counter,": ",replicationNumber,"/",replications," ",populationModelIndex,"/27 ",sampleSizeIndex,"/3",sep=""))
			
			# Which rows of the desing matrix are we running now
			startIndex<-populationModelIndex*27+sampleSizeIndex*9-35
			endIndex<-startIndex+8

			
			# Write the population model, tested model and 
			# simulation paramaters to a file that will be used as 
			# input for MapReduce
			
			# The protocol is to have data as a single line as tab
			# separated in the following order:
			# Counter, Index of first row, Index of last row
			
			cat(counter,startIndex,endIndex,sep="\t",file="input.txt",append=TRUE)

			cat("\t",file="input.txt",append=TRUE)

			# Then the population model 

			cat( populationModel$covariances,populationModel$paths,sep="\t",file="input.txt",append=TRUE)

			# Then the three tested models
			
			for(omittedPathsShareIndex in 1:3){
				
				extraPathsIndex <- designMatrix[startIndex+omittedPathsShareIndex*3-3,6]
	
				cat("\t",file="input.txt",append=TRUE)
				cat(testedModels[[(omittedPathsShareIndex-1)*3+extraPathsIndex]],sep="\t",file="input.txt",append=TRUE)

			}
			
			cat("\n",file="input.txt",append=TRUE)
		}
	}
}
