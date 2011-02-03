#! /usr/bin/env Rscript

# This file contains all function definitions needed for combining the results from MapReduce


#
# Reads a matrix or a data frame from text connection. This function is needed
# because calling read.delim directly will read to the end of the connection
# instead of stopping to first blank line.
#

readData <- function(connection,allow.na=FALSE){
	
	# Create a file and write to it until we hit a blank line
	
	file<-file("")
	
	while (length(line <- trim(readLinesWithCounter(connection, n = 1, warn = FALSE))) > 0) {
		debugPrint(line)
		if(line==""){
			writeLines(line,con=file)
			break()
		}
		else if(grepl("^\\[[0-9]+\\]",line)){
			# Just skip any lines created with Print
		}
		else{
			writeLines(line,con=file)
		}
	}
	
	res<-read.delim(file)
	if(! allow.na){
		if(sum(is.na(res))>0){
			print(res)
			stop(paste("Matrix contains NA values. Line number",lineNumber))
		}
	}
	return(res)
}

readLinesWithCounter <- function(con, n = 1, warn = FALSE){
	lineNumber<<-lineNumber+1
	debugPrint(lineNumber)
	return(readLines(con, n , warn ))
}