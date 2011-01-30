#! /usr/bin/env Rscript
#
# This file reads the results from the MapReduce from stdin, parses the results
# and saves these as an R workspace image that can be later loaded for analyses.
#

source("parameters.R")
source("functions.R")

# Initialize matrices
designMatrix<-createdesignMatrix()

# Paths are stored as a matrix or data.frame (which ever read.table defaults to) 

paths<-NULL

# Three dimensional array. First is the replication number, second is the
# design number and third is the estimation methdod

correlations<-array(list(rep(NULL,729*replications*4)),dim=c(replications,729,4),dimnames=c("replication","designNumber","analysis"))

populationFactorLoadings<-array(list(rep(NULL,729*replications)),dim=c(replications,729),dimnames=c("replication","designNumber",))

constructEstimateSdsByModel<-array(list(rep(NULL,729*replications*4)),dim=c(replications,729,4),dimnames=c("replication","designNumber","analysis"))

constructEstimateSdsByData<-array(list(rep(NULL,729*replications*4)),dim=c(replications,729,4),dimnames=c("replication","designNumber","analysis"))


con <- file("stdin", open = "r")

while (length(line <- readLines(con, n = 1, warn = FALSE)) > 0) {
	
	# Detect running a new design
	if(line=="Design"){
	
		line <- readLines(con, n = 1, warn = FALSE)

		specification<-sapply(unlist(strsplit(line, split="\t")),as.numeric)
		
		# Loop over the design numbers and store the results to the results 
		# objects
			
		constructCount<-sqrt((length(specification)-3)/2)
		populationPaths<-matrix(specification[c(4+constructCount^2:length(specification))],ncol=constructCount)
			
		populationPathData<-NULL
			
		names<-paste("C",c(1:constructCount),sep="")
			
		for(i in 1:constructCount){
			for(j in 1:constructCount){
				if( ! is.na(populationPaths[i,j]) & populationPaths[i,j]!=0){
					# "From","To","Estimate","Mean.Boot","Std.Error","perc.05","perc.95"
					populationPathData<-rbind(populationPathData,c(names[j],names[i],populationPaths[i,j],NA,NA,NA,NA)
				}
			}
		}

		for(designNumber in specification[2]:specification[3]){

			#Add population paths
			paths<-rbind(paths,cbind(populationPathData,replication,designNumber,"truevalue"))

		}

		#Reset all data objects

		analysis<-NULL
	}

	# detect that we will read analysis results next
	
	else if (line=="sumscale" | line=="component" | line=="factor" | line=="pls"){
		# If analysis is null, this means that we just finished reading
		# population data
		
		
		analysis<-line
	}
	
	# If analysis is not specified yet, we are reading data related to the sample
	
	else if (is.null(analysis)){

		# Detect that the next thing is results for different data
		if(grep("Data [123]",line)){
			dataNumber<-as.numeric(substr(line,5,6))
		}
		else if(line=="Population factor loadings"){
			tempLoadings[[dataNumberl]]<-read.table("stdin",sep="\t")

			#Add factorloadings
			
			for(designNumber in specification[2]:specification[3]){
				
				# 7th column of the desing matrix is the identity column for 
				# data. Check that it matches and record the factor loadings
				
				if(designMatrix[desingNumber,7]==dataNumber){
					populationFactorLoadings[[replication,designNumber]]<-tempLoadings
				}
			}
		}
	}

	# Read analysis results

	else {
		# Read measurement indeterminancy things
		if(grep("Data [123]",line)){
			dataNumber<-as.numeric(substr(line,6,6))
			sds<-read.table("stdin",sep="\t")
			
			for(designNumber in specification[2]:specification[3]){
				
				# 7th column of the desing matrix is the identity column for 
				# data. Check that it matches and record the standard deviations
				
				if(designMatrix[desingNumber,7]==dataNumber){
					constructEstimateSdsByData[[replication,designNumber]]<-sds
				}
			}

			
		}
		else if(grep("Model [123]",line)){
			modelNumber<-as.numeric(substr(line,6,6))
			sds<-read.table("stdin",sep="\t")
			
			for(designNumber in specification[2]:specification[3]){
				
				# 5th column of the desing matrix is the identity column for 
				# model. Check that it matches and record the standard 
				# deviations
				
				if(designMatrix[desingNumber,5]==modelNumber){
					constructEstimateSdsByModel[[replication,designNumber]]<-sds
				}
			}
		}
		# Read results from a replication
		else if(grep("Row index [123]+",line)){
			designNumber<-as.numeric(substr(line,11,nchar(line)))
		}
		else if(line=="Paths"){
			newPathSet<-read.table("stdin",sep="\t")
			# Add information about the analys and append to the results
			paths<-rbind(paths,cbind(newPathSet,replication,designNumber,analysis)
		}
		else if(line=="Correlations"){
			correlations[[replication,designNumber,analysis]]<-read.table("stdin",sep="\t")
		}

	}
}

# Save the data that were parsed from stdin

save(c("designMatrix","paths","correlations","populationFactorLoadings","constructEstimateSdsByModel","constructEstimateSdsByData"))
,file="data.RData")