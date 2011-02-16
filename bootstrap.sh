#!/bin/bash
#
# This file installs the packages needed for the MapReduce job.
#


# debian R upgrade (Source: http://www.r-bloggers.com/bootstrapping-the-latest-r-into-amazon-elastic-map-reduce/)
# Ensure that we are only using stable packages by over writing the apt sources and apt preferences

echo "deb http://ftp.us.debian.org/debian stable main non-free contrib" | sudo tee /etc/apt/sources.list

echo 'Package: *
Pin: release a=stable
Pin-Priority: 990' | sudo tee /etc/apt/preferences

cat /etc/apt/sources.list
cat /etc/apt/preferences

sudo apt-get -t stable update 
sudo DEBIAN_FRONTEND=noninteractive apt-get -t stable install --yes --force-yes r-recommended


# plspm: The PLS package
# psych: Needed to run regressions using correlation matrix as data
# rest: dependencies for plspm

# The dependencies are installed manually because they are really not needed and we do not want
# to pull the entire dependency tree behind these packages.

# Running with sudo to make sure that R has permissions to install packages

echo 'install.packages(c("plspm","amap","diagram","shape","psych"),repos="http://cran.r-project.org")' | sudo R --no-save 

