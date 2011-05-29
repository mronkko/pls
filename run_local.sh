#!/bin/bash

#
# This file runs the analysis scripts on local computer
#

PATH=./R-2.13.0/bin:$PATH

#Check that the R is correct
which R
R --version

THREADS=20
COUNTER=0

while read line; do 
	let COUNTER=COUNTER+1
	echo $COUNTER
	#Check that there are less threads running than the maximum limit
	

	while [ `ps -C R | wc -l` -gt $THREADS ]; do
		echo "More than $THREADS R processes running. Wait for some of them to exit before proceeding"
		sleep 5
	done

	#Check that an output file does not already exist
	
	if [ ! -e output/$COUNTER.csv ]; then
		echo $line | nice ./reduce.R > output/$COUNTER.csv &
	fi
done < input.txt