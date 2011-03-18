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


elastic-mapreduce --create --stream --input s3n://$INPUT/input.txt --output s3n://$OUTPUT/  --mapper s3n://$INPUT/map.R --reducer s3n://$INPUT/reduce.R  --cache-archive s3n://$INPUT/cache.zip#include --log-uri s3://$LOGS/ --instance-type m1.small --instance-count c1.medium  --bootstrap-action s3n://$INPUT/bootstrap.sh  --enable-debugging --ssh --bootstrap-action s3://elasticmapreduce/bootstrap-actions/configure-hadoop --args "-m,mapred.reduce.tasks=`cat input.txt | wc | sed 's/ //'`"


#
# Instance types that can be used (these provide the most bang for buck)
# 
# m1.small
# c1.medium
# c1.xlarge
#
# These contain more memory 
#
# m1.large
#
# Path to log files on Hadoop master node is
# /mnt/var/lib/hadoop/mapred/taskTracker/archive/mronkko-pls-analysisfiles/
#