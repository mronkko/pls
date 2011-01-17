#! /usr/bin/env Rscript

#
# Prepares an input file for the MapReduce jobs. The script runs several
# replicatiosn and ror each replication, it generates
# one population model for each combination of population model parametrization
# (number of constructs, paths set to zero, population path values). The it
# generates all combinations of tested model parameters (paths omitted, extra
# paths) and passes the population model and tested model combinations to
# the reducer that will create datasetss based on the population model and
# evaluate the tested model with these data.
#
#


# Read simulation parameters
source("parameters.R")

# Read function definitions
source("functions.R")

# Keep a counters of reduce jobs and model numbers
counter<-0
populationModelNumber<-0
testedModelNumber<-0

for(replicationNumber in 1:replications){
    # Generate all population models and all tested models in mapper and send these to reducer that will generate the construct and indicator scores and test the models.
    
    for(numberOfConstructsIndex in 1:length(numberOfConstructs)){
	    for(expectedNumberOfOutgoingPathsIndex in 1:length(expectedNumberOfOutgoingPaths)){
		    for(populationPathValuesIndex in 1:length(populationPathValues)){
		    
		    	# Generate a population model here. We need two matrices. 
		    	# First matrix contains 1s and 0s to store which pats are specified
		    	# Second matrix contains the path values
		    
		    	populationModelWhichPaths <- generateRandomModel(numberOfConstructs[numberOfConstructsIndex],expectedNumberOfOutgoingPaths[expectedNumberOfOutgoingPathsIndex])

				populationModel<-setPopulationModelPathValues(populationModelWhichPaths,populationPathValues[[populationPathValuesIndex]])
				
				populationModelNumber <- populationModelNumber + 1
		    
			    for(omittedPathsShareIndex in 1:length(omittedPathsShare)){
				    for(extraPathsIndex in 1:length(extraPaths)){
				    
					    testedModel <- generateTestedModel(populationModelWhichPaths,omittedPathsShare[omittedPathsShareIndex],extraPaths[extraPathsIndex])

#						print("population")
#						print(populationModel)
#						print("tested")
#						print(testedModel)
						
						testedModelNumber <- testedModelNumber + 1

						counter <- counter +1

						#write a progress indicator on screen

						print(paste(counter,"-",replicationNumber,"/",replications," ",numberOfConstructsIndex,"/",length(numberOfConstructs)," ",expectedNumberOfOutgoingPathsIndex,"/",length(expectedNumberOfOutgoingPaths)," ",populationPathValuesIndex,"/",length(populationPathValues)," ",omittedPathsShareIndex,"/",length(omittedPathsShare)," ",extraPathsIndex,"/",length(extraPaths),sep=""))


						# Write the population model, tested model and 
						# simulation paramaters to a file that will be used as 
						# input for MapReduce
						
						# The protocol is to have data as a single line as tab
						# separated in the following order:
						# Counter
						# Population model number
						# Tested model number
						# 5 variables that describe the population model and tested model parameters
						# Population model
						# Tested model

						cat(counter, populationModelNumber, testedModelNumber,
						numberOfConstructs[numberOfConstructsIndex],expectedNumberOfOutgoingPaths[expectedNumberOfOutgoingPathsIndex],unlist(populationPathValues[populationPathValuesIndex]),
						omittedPathsShare[omittedPathsShareIndex],extraPaths[extraPathsIndex], populationModel$covariances,
						populationModel$paths,testedModel, sep="\t",file="input.txt",append=TRUE)
						cat("\n",file="input.txt",append=TRUE)
					}
				}
			}
		}
	}
}