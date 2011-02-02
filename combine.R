#! /usr/bin/env Rscript
#
# This file reads the results from the MapReduce from stdin, parses the results
# and saves these as an R workspace image that can be later loaded for analyses.
#

source("include/parameters.R")
source("include/functions.R")
source("include/functionsCombine.R")

# Needed for trim-function
library(R.oo) 

# Initialize matrices
designMatrix<-createdesignMatrix()

# Paths are stored as a data.frame (which ever read.table defaults to) 

paths<-NULL

# Three dimensional array. First is the replication number, second is the
# design number and third is the estimation methdod

correlations<-array(list(rep(NULL,729*replications*4)),dim=c(replications,729,4),dimnames=c("replication","designNumber","analysis"))

populationFactorLoadings<-array(list(rep(NULL,729*replications)),dim=c(replications,729),dimnames=c("replication","designNumber"))

models<-array(list(rep(NULL,729*replications*2)),dim=c(replications,729,2),dimnames=c("replication","designNumber","type"))


constructEstimateSdsByModel<-array(list(rep(NULL,729*replications*4)),dim=c(replications,729,4),dimnames=c("replication","designNumber","analysis"))

constructEstimateSdsByData<-array(list(rep(NULL,729*replications*4)),dim=c(replications,729,4),dimnames=c("replication","designNumber","analysis"))


con <- file("results.txt", open = "r")


#Reset the analysis object
analysis<-NULL

while (length(line <- trim(readLines(con, n = 1, warn = FALSE))) > 0) {
	
	debugPrint(line)
	
	#If the line is empty, do nothing
	if (line==""){
		next()
	}
	# Detect running a new design
	else if(line=="Design"){
	
		line <- readLines(con, n = 1, warn = FALSE)

		specification<-sapply(unlist(strsplit(line, split="\t")),as.numeric)
		replication<-specification[1]
		
		print(paste("Replication ",specification[1],", designs  ",specification[2],"-",specification[3],sep=""))
		
		# Loop over the design numbers and store the results to the results 
		# objects
			
		constructCount<-sqrt((length(specification)-3)/5)
		
		populationPaths<-matrix(specification[c((4+constructCount^2):(3+2*constructCount^2))],ncol=constructCount)

		populationPathData<-NULL
			
		names<-paste("C",c(1:constructCount),sep="")
			
		for(i in 1:constructCount){
			for(j in 1:constructCount){
				if( ! is.na(populationPaths[i,j]) & populationPaths[i,j]!=0){
					# "From","To","Estimate","Mean.Boot","Std.Error","perc.05","perc.95"
					
					populationPathData<-rbind(populationPathData,data.frame(names[j],names[i],populationPaths[i,j],NA,NA,NA,NA))

					
				}
			}
		}

		if(! is.null(populationPathData)){
			for(designNumber in specification[2]:specification[3]){

				#Add population paths
				newRow<-cbind(populationPathData,replication,designNumber,"truevalue")
				colnames(newRow)<-c("From","To","Estimate","Mean.Boot","Std.Error","perc.05","perc.95","replication","designNumber","analysis")
				
				paths<-rbind(paths,as.data.frame(newRow))
				
			}
		}
		
		
		
		# Store population model and tested models
		
		for(designNumber in specification[2]:specification[3]){
			models[[replication,designNumber,1]]<- ! (is.na(populationPaths) | populationPaths==0)

			models[[replication,designNumber,2]]<- matrix(specification[c((4+constructCount^2*(designMatrix[designNumber,5]+1)):(3+constructCount^2*(designMatrix[designNumber,5]+2)))],ncol=constructCount)
		}

		# Set the column names. These are not set
		
		analysis<-NULL
	}

	# detect that we will read analysis results next
	
	else if (line %in% analysisTypes){
		# If analysis is null, this means that we just finished reading
		# population data
		
		
		analysis<-which(analysisTypes==line)
	}
	
	# If analysis is not specified yet, we are reading data related to the sample
	
	else if (is.null(analysis)){

		# Detect that the next thing is results for different data

		if(grepl("Data [123]",line)){
			dataNumber<-as.numeric(substr(line,5,6))
		}
		else if(line=="Population factor loadings"){

			
			tempLoadings<-readData(con)

			#Add factorloadings
			for(designNumber in specification[2]:specification[3]){
				
				# 7th column of the design matrix is the identity column for 
				# data. Check that it matches and record the factor loadings
				
				if(designMatrix[designNumber,7]==dataNumber){
					populationFactorLoadings[[replication,designNumber]]<-tempLoadings
				}
			}
		}
	}

	# Read analysis results

	else {
		# Read measurement indeterminancy things
		if(grepl("Data [123]",line)){
			dataNumber<-as.numeric(substr(line,6,6))
			sds<-readData(con)
			
			for(designNumber in specification[2]:specification[3]){
				
				# 7th column of the design matrix is the identity column for 
				# data. Check that it matches and record the standard deviations
				
				if(designMatrix[designNumber,7]==dataNumber){
					constructEstimateSdsByData[[replication,designNumber,analysis]]<-sds
				}
			}

			
		}
		else if(grepl("Model [123]",line)){
			modelNumber<-as.numeric(substr(line,7,7))
			sds<-readData(con)
			
			for(designNumber in specification[2]:specification[3]){
				
				# 5th column of the design matrix is the identity column for 
				# model. Check that it matches and record the standard 
				# deviations
				
				if(designMatrix[designNumber,5]==modelNumber){
					constructEstimateSdsByModel[[replication,designNumber,analysis]]<-sds
				}
			}
		}
		# Read results from a replication
		else if(grepl("Row index [0-9]+",line)){
			designNumber<-as.numeric(substr(line,11,nchar(line)))
		}
		else if(line=="Paths"){
			newPathSet<-readData(con)
			# Add information about the analys and append to the results
			
			newRow<-cbind(newPathSet,replication,designNumber,analysisTypes[analysis])
			colnames(newRow)<-c("From","To","Estimate","Mean.Boot","Std.Error","perc.05","perc.95","replication","designNumber","analysis")
			

			paths<-rbind(paths,newRow)


		}
		else if(line=="Correlations"){
			correlations[[replication,designNumber,analysis]]<-readData(con)
		}

	}
}

# Save the data that were parsed from stdin

save(designMatrix,paths,correlations,populationFactorLoadings,constructEstimateSdsByModel,constructEstimateSdsByData,models
,file="data.RData")