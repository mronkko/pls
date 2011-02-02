#! /usr/bin/env Rscript
#
# This file installs the packages needed for the MapReduce job.
#

# plspm: The PLS package
# MASS: Needed for generating data from covariance matrices
# R.oo: Needed for trim-function
install.packages(c("plspm","MASS","R.oo"),repos="http://cran.r-project.org",lib=".")



