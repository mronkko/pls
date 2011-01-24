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


# Read simulation parameters
source("parameters.R")

# Read function definitions
source("prepareFunctions.R")


#
# We start by generating a fractional factorial desing matrix.
#
# See the this article for details:
# Xu, H. 2005. “A catalogue of three-level regular fractional factorial designs,” Metrika (62:2), pp. 259–281.
#
# The design matrix is 11-6.1 from the 234 replications set from the article
# above.
#

generator <- matrix(c(	1,0,0,1,
						0,1,1,2),ncol=4,byrow=TRUE)

for(i in 3:5){
	generator<-cbind(rbind(generator,0),c(rep(0,i-1),1),rbind(generator,1),rbind(generator,2))
}

#
# permutations of a b c d e 
#
permutations<-matrix(c(c(0:242)%%3,c(0:242)%/%3%%3,c(0:242)%/%9%%3,c(0:242)%/%27%%3,c(0:242)%/%81),ncol=5)

# Design matrix is the product of permutations and generator matrix, mod 3

fullDesignMatrix<- ( permutations %*% generator ) %% 3

# Get the indices from the 11-6.1 design reported in Xu, 2005

print(ncol(fullDesignMatrix))
print(nrow(fullDesignMatrix))

desingMatrix<-fullDesignMatrix[,c(1,2,5,14,41,63,27,72,79,93,114)]

debugPrint(desingMatrix)

stop("TODO: Implement the design matrix thing here")

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
	    
	    
	    	# Generate all tested models for this population model
	    	
	    	testedModels<-NULL
	    	
		    for(omittedPathsShareIndex in 1:length(omittedPathsShare)){
			    for(extraPathsIndex in 1:length(extraPaths)){
			    
				    testedModels <-c(testedModels, generateTestedModel(populationModelWhichPaths,omittedPathsShare[omittedPathsShareIndex],extraPaths[extraPathsIndex]))
				}
			}
			
			# Loop over data related conditions
			
				for(sampleSizeIndex in 1:length(sampleSizes)){
					for(indicatorCountIndexs in 1:length(indicatorCounts)){
						for(factorLoadingIndex in 1:length(factorLoadings)){
							for(factorLoadingIntervalIndex in 1:length(factorLoadingIntervals)){
								for(maxErrorCorrelationIndex in 1:length(maxErrorCorrelations)){
									for(methodVarianceIndex in 1:length(methodVariances)){

										counter <- counter +1

										#write a progress indicator on screen
										print(paste(counter,"-",replicationNumber,"/",replications," ",numberOfConstructsIndex,"/",length(numberOfConstructs)," ",expectedNumberOfOutgoingPathsIndex,"/",length(expectedNumberOfOutgoingPaths)," ",populationPathValuesIndex,"/",length(populationPathValues)," ",sampleSizeIndex,"/",length(sampleSizes)," ",indicatorCountIndexs,"/",length(indicatorCounts)," ",factorLoadingIndex,"/",length(factorLoadings)," ",factorLoadingIntervalIndex,"/",length(factorLoadingIntervals)," ",maxErrorCorrelationIndex,"/",length(maxErrorCorrelations)," ",methodVarianceIndex,"/",length(methodVariances),sep=""))
										
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

										cat(counter, populationModelNumber, numberOfConstructsIndex,expectedNumberOfOutgoingPathsIndex,populationPathValuesIndex,sampleSizeIndex, indicatorCountIndexs, factorLoadingIndex, factorLoadingIntervalIndex, maxErrorCorrelationIndex, methodVarianceIndex, populationModel$covariances,populationModel$paths,testedModels, sep="\t",file="input.txt",append=TRUE)

										cat("\n",file="input.txt",append=TRUE)
									}
								}
							}
						}
					}
				}
			}
		}
	}
}
