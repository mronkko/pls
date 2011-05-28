#!/bin/bash

#
# This file runs the analysis scripts on local computer
#

THREADS=20
COUNTER=0

for line in $(cat input.txt); do 
	let COUNTER=COUNTER+1
	
	#Check that there are less threads running than the maximum limit
	

	while [ `ps -C R | wc -l` -gt $THREADS ]; do
		echo "More than $THREADS R processes running. Wait for some of them to exit before proceeding"
		sleep 5
	done	

	#Check that an output file does not already exist
	
	if [ ! -e output/$COUNTER.csv ]; do
		echo $line
		echo $line | ./reduce.R > output/$COUNTER.csv &
	done
done
