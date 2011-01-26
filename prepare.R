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
# We start by generating a fractional factorial desing matrix. Running this 
# simulation with full factorial desing is not feasible due to the number of 
# factors.
#
# See the this article for details:
# Xu, H. 2005. “A catalogue of three-level regular fractional factorial
# designs,” Metrika (62:2), pp. 259–281.
#
# The key problem with creating the design is that the simulation runs are 
# nested: We have several sets of construct scores. We want to run at least
# three different measurements for each score while keeping the tested model
# constant to estimate test-retest stability of the estimated (Constraint 1).
# Also we need to use each indicator data to estimate at least three different 
# models so that we can analyze the stability of the measurement (Constraint 2).
#
# All replications must be independent, so these two constraints need to be
# taken care of in the design matrix. 
#
# We should have one identity column that does not affect the construt scores to # be able to estimate three models for each set of constucts. This identity 
# column should not affect the indicator columns either to keep both the
# indicator data constant (no changes in constructs or measurements) over three # model estimations
#
# The first desing that fills this criteria is 11-5.10 from the 729 replications # set from the article above.
#
# See http://www.stat.ucla.edu/~hqxu/pub/ffd/ffd729.pdf
#


#
# 729 runs would allow 14 factors with resolution V (14-8.1)

generator <- matrix(c(	1,0,0,1,
						0,1,1,2),ncol=4,byrow=TRUE)

for(i in 3:6){
	generator<-cbind(rbind(generator,0),c(rep(0,i-1),1),rbind(generator,1),rbind(generator,2))
	
}


#
# This command prints the desing matrix. The output is included below the 
# command.
#
# print(generator[,c(1,2,5,14,41,122,63,149,166,188,222)])
#
#      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11]
# [1,]    1    0    0    0    0    0    0    1    0     0     0
# [2,]    0    1    0    0    0    0    1    2    1     1     0
# [3,]    0    0    1    0    0    0    1    2    0     2     1
# [4,]    0    0    0    1    0    0    1    1    0     1     1
# [5,]    0    0    0    0    1    0    1    0    1     1     2
# [6,]    0    0    0    0    0    1    0    1    1     1     1
#
# The factors are assigned as follows
# Tested model: 1, 8
# Measurement: 2,7,9,10,11
# Population model and true scores: 3,4,5,6
#

#
# permutations of a b c d e f
#
permutations<-matrix(c(c(0:728)%%3,c(0:728)%/%3%%3,c(0:728)%/%9%%3,c(0:728)%/%27%%3,c(0:728)%/%81%%3,c(0:728)%/%243),ncol=6)

print(permutations)

# Design matrix is the product of permutations and generator matrix, mod 3

fullDesignMatrix<- ( permutations %*% generator ) %% 3

# Get the indices from the 11-5.10 design reported in the online appendix of Xu,
# 2005

desingMatrix<-fullDesignMatrix[,c(1,2,5,14,41,122,63,149,166,188,222)]


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
