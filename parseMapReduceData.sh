#!/bin/bash
#
# Reads in the part-xxxxx files produced by the MapReduce algorithm and combines these into
# two CSV files that
#

# Create csv headeds headers

echo '"type"	"replication"	"designNumber"	"analysis"	"to"	"from"	"trueCorrelation"	"estimatedCorrelation"	"correlationAttenuationCoefficient"	"regressionTrueScore"	"regressionEstimate"	"regressionSE"' > data/relationships.csv
echo '"type"	"replication"	"designNumber"	"analysis"	"construct"	"CR"	"AVE"	"minFactorLoading"	"meanFactorLoading"	"maxCrossLoading"	"maxCorrelationWithOtherConstruct"	"trueScoreCorrelation"	"deltaR2constructs"	"deltaR2errors"	"estimatedR2"	"trueR2"	"trueConstructR2"	"estimatedConstructR2"	"sdByData"	"sdByModels"	"incomingPathsCorrect"	"incomingPathsExtra"	"incomingPathsOmitted"	"outgoingPathsCorrect"	"outgoingPathsExtra"	"outgoingPathsOmitted"' > data/constructs.csv
echo '"type"	"replication"	"designNumber"	"analysis"	"SRMR"	"GlobalGoF"	"meanSquareResiduals"	"varianceExplainedCommonFactor"	"smallestPositiveCorrelation"	"secondSmallestPositiveCorrelation"' > data/models.csv

# Loop over the part-xxxxx files and print the content to the approriate filrs

find output  -exec grep '^M' {} \; >> data/models.csv
find output  -exec grep '^C' {} \; >> data/constructs.csv
find output  -exec grep '^R' {} \; >> data/relationships.csv

