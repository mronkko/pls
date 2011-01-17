#! /usr/bin/env Rscript
#
# This file is the mapper of MapReduce.
#


# Read the stdin for which replications to run
con <- file("stdin", open = "r")


while (length(line <- readLines(con, n = 1, warn = FALSE)) > 0) {
    cat(line)
}
close(con)