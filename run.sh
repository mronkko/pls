#!/bin/bash

#
# These are the buckets that are used. They must exist before the script is used
#

LOGS=mronkko-pls-logs
INPUT=mronkko-pls-analysisfiles
OUTPUT=mronkko-pls-results

#
# Copy the files to S3
#

zip -j cache.zip include/functions.R include/parameters.R

s3cmd put cache.zip input.txt map.R reduce.R bootstrap.sh s3://$INPUT/


#
# Sets up a map reduce job on Amazon cloud
#


elastic-mapreduce --create --stream --input s3://$INPUT/input.txt --output s3://$OUTPUT/  --mapper s3://$INPUT/map.R --reducer s3://$INPUT/reduce.R  --cache-archive s3://$INPUT/cache.zip#include --log-uri s3://$LOGS/ --instance-type c1.medium --instance-count 2  --bootstrap-action s3://$INPUT/bootstrap.sh


#
# instance types that can be used (these provide the most bang for buck)
# 
# m1.small
# c1.medium
# c1.xlarge
#
# 

#
# Path to log files on Hadoop master node is
# /mnt/var/lib/hadoop/mapred/taskTracker/archive/mronkko-pls-analysisfiles/
#