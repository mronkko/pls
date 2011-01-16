#! /usr/bin/env Rscript
# The reducer receives a population model and tested model. It runs PLS and summed scales on all specified experimental conditions for these models and outputs the results as CSV. 

source("parameters.R")

options(warn=-1)

con <- file("stdin", open = "r")

#One line means one population model-tested model pair. See map.R for details on the protocol

while (length(line <- readLines(con, n = 1, warn = FALSE)) > 0) {

}

close(con)
