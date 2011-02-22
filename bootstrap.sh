#!/bin/bash
#
# This file installs the packages needed for the MapReduce job.
#


# debian R upgrade (Source: http://www.r-bloggers.com/bootstrapping-the-latest-r-into-amazon-elastic-map-reduce/)

# Ensure that we are only using stable packages by over writing the apt sources and apt preferences. Add CRAN repository that contains R 2.11 backported for Debian Lenny that Amazon uses.

#echo 'deb http://http.us.debian.org/debian   lenny         main contrib non-free
#deb http://security.debian.org         lenny/updates main contrib non-free
#deb http://ftp.heanet.ie/mirrors/cran.r-project.org/bin/linux/debian lenny-cran/' | sudo tee /etc/apt/sources.list

#sudo apt-get update 
#sudo DEBIAN_FRONTEND=noninteractive apt-get install --yes --force-yes r-recommended

# Download and install R version 2.12.1 that has been compiled for this computer

wget http://www.tkk.fi/u/mronkko/R.zip
unzip R.zip 
cd R-2.12.1/
sudo make install

# Check that installing succceeded
which R
R --version 

# plspm: The PLS package
# psych: Needed to run regressions using correlation matrix as data
# rest: dependencies for plspm

# The dependencies are installed manually because they are really not needed and we do not want
# to pull the entire dependency tree behind these packages.

# Running with sudo to make sure that R has permissions to install packages

echo 'install.packages(c("plspm","amap","diagram","shape","psych"),repos="http://cran.r-project.org")' | sudo R --no-save 

