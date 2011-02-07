#!/bin/bash
#
# Reads in the part-xxxxx files produced by the MapReduce algorithm and combines these into
# two CSV files that
#

# Create csv headeds headers

echo '"type"	"replication"	"designNumber"	"analysis"	"to"	"from"	"trueCorrelation"	"estimatedCorrelation"	"correlationAttenuationCoefficient"	"regressionTrueScore"	"regressionEstimate"	"regressionSE"' > data/relationships.csv
echo '"type"	"replication"	"designNumber"	"analysis"	"construct"	"CR"	"AVE"	"minFactorLoading"	"meanFactorLoading"	"maxCrossLoading"	"maxCorrelationWithOtherConstruct"	"trueScoreCorrelation"	"deltaR2"	"estimatedR2"	"trueR2"	"sdByData"	"sdByModels"	"incomingPathsCorrect"	"incomingPathsExtra"	"incomingPathsOmitted"	"outgoingPathsCorrect"	"outgoingPathsExtra"	"outgoingPathsOmitted"' > data/constructs.csv

# Loop over the part-xxxxx files and print the content to the approriate filrs

find data -name 'part-*' -exec grep '^C' {} \; >> data/constructs.csv
find data -name 'part-*' -exec grep '^R' {} \; >> data/relationships.csv

