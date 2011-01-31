#! /usr/bin/env Rscript

# This file contains all function definitions needed for combining the results from MapReduce


#
# Reads a matrix or a data frame from text connection. This function is needed
# because calling read.delim directly will read to the end of the connection
# instead of stopping to first blank line.
#

readData <- function(connection){
	
	# Create a fifo and write to it until we hit a blank line
	
	fifo<-fifo("")
	
	debugPrint("reading data....")
	
	while (length(line <- readLines(connection, n = 1, warn = FALSE)) > 0) {
		if(line==""){
			break()
		}
		else if(grepl("^\\[[0-9]+\\]",line)){
			# Just skip any lines created with Print
		}
		else{
			writeLines(line,con=fifo)
		}
	}
	
	res<-read.delim(fifo)
	debugPrint(res)
	return(res)
}