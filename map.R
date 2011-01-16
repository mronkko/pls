#! /usr/bin/env Rscript
#
# This file is the mapper of MapReduce. It reads in replication numbers from
# stdin and does these replications. For each replication, the mapper generates
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

# Read the stdin for which replications to run
con <- file("stdin", open = "r")

# Keep a counters of reduce jobs and model numbers
counter<-0
populationModelNumber<-0
testedModelNumber<-0

while (length(line <- readLines(con, n = 1, warn = FALSE)) > 0) {
    
    # Generate all population models and all tested models in mapper and send these to reducer that will generate the construct and indicator scores and test the models.
    
    for(numberOfConstructsIndex in 1:length(numberOfConstructs)){
	    for(expectedNumberOfOutgoingPathsIndex in 1:length(expectedNumberOfOutgoingPaths)){
		    for(populationPathValuesIndex in 1:length(populationPathValues)){
		    
		    	# Generate a population model here. We need two matrices. 
		    	# First matrix contains 1s and 0s to store which pats are specified
		    	# Second matrix contains the path values
		    
		    	populationModelWhichPaths <- generateRandomModel(numberOfConstructs[numberOfConstructsIndex],expectedNumberOfOutgoingPaths[expectedNumberOfOutgoingPathsIndex])

				populationModel<-setPopulationModelPathValues(populationModelWhichPaths,,populationPathValues[populationPathValuesIndex])
				
				populationModelNumber <- populationModelNumber + 1
		    
			    for(omittedPathsShareIndex in 1:length(omittedPathsShare)){
				    for(extraPathsIndex in 1:length(extraPaths)){
				    
					    testedModel <- generateTestedModel(populationModelWhichPaths,omittedPathsShare[omittedPathsShareIndex],extraPaths[extraPathsIndex])

						testedModelNumber <- testedModelNumber + 1

						counter <- counter +1

						# Pass the population model and the
						# The protocol is to have data as a single line as tab separated in the following order:
						# Counter
						# Population model number
						# Tested model number
						# 5 variables that describe the population model and tested model parameters
						# Population model
						# Tested model

						cat(paste(counter, populationModelNumber, testedModelNumber,numberOfConstructs[numberOfConstructsIndex],expectedNumberOfOutgoingPaths[expectedNumberOfOutgoingPathsIndex],populationPathValues[populationPathValuesIndex],omittedPathsShare[omittedPathsShareIndex],extraPaths[extraPathsIndex], populationModel, testedModel, sep="\t"),"\n", sep="") 

					}
				}
			}
		}
	}

}
close(con)