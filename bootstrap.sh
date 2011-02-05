#!/bin/bash
#
# This file installs the packages needed for the MapReduce job.
#

uname -a

# plspm: The PLS package
# psych: Needed to run regressions using correlation matrix as data
# rest: dependencies for plspm

# The dependencies are installed manually because they are really not needed and we do not want
# to pull the entire dependency tree behind these packages.

# Running with sudo to make sure that R has permissions to install packages

echo 'install.packages(c("plspm","amap","diagram","shape","psych"),repos="http://cran.r-project.org")' | sudo R --no-save 

