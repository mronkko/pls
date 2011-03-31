#!/bin/bash

#
# These are the S3 buckets that are used. They must exist before the script is 
# used. Change these to buckets that you own before running this file.
#

LOGS=logs
INPUT=analysisfiles
OUTPUT=results

#
# Copy the files to S3
#

zip -j cache.zip include/functions.R include/parameters.R include/plspm.R include/plspm-internal.R

# You must have s3cmd properly installed.

s3cmd put cache.zip input.txt map.R reduce.R bootstrap.sh s3://$INPUT/


#
# Sets up a map reduce job on Amazon cloud. The embedded shell script sets the 
# number of reduce tasks to equal the number of input file. This is not 
# nescessary, but it will speed up getting the first results.
#


elastic-mapreduce --create --stream --input s3n://$INPUT/input.txt --output s3n://$OUTPUT/  --mapper s3n://$INPUT/map.R --reducer s3n://$INPUT/reduce.R  --cache-archive s3n://$INPUT/cache.zip#include --log-uri s3://$LOGS/ --instance-type c1.medium --instance-count 20  --bootstrap-action s3n://$INPUT/bootstrap.sh  --bootstrap-action s3://elasticmapreduce/bootstrap-actions/configure-hadoop --args "-m,mapred.reduce.tasks=`cat  input.txt | wc -l | sed 's/^[^0-9]*//'`"


#
# Instance types that can be used (these provide the most bang for the buck)
# 
# m1.small
# c1.medium
# c1.xlarge
#