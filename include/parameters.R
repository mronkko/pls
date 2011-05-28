# Overall study related parameters

replications<-500

# Model related parameters ( 3 x 3 x 3 = 27 combinations)
# 3: 1-3
numberOfConstructs <-c(4,8,12)
expectedNumberOfOutgoingPaths <-c(1,2,3)
populationPathValues<-list(c(-.5,.5),c(0,.5),c(0,0))

# Parameters related to tested models (3 x 3 = 9 combinations)
# 2: 4-5
omittedPathsShare<-c(0,.25,.5)
extraPaths<-c(0,1,2)


# Data related parameters (3 combinations)
# 1: 6
sampleSizes<-c(50,100,250)

# Measurement model related parameters. These all are population values. ( 3 x 3 x 3 x 3 x 3 = 243 combinations)
# 5: 7-11
indicatorCounts<-c(3,4,6)
factorLoadings<-c(.4,.6,.8)
factorLoadingIntervals<-c(0,.1,.2)
maxErrorCorrelations<-c(0,.2,.4)
methodVariances<-c(0,.1,.2)

# Do not change this. 
analysisTypes<-c("sumscale","component","factor","pls_Standard","pls_IndividualSignChanges","pls_ConstructLevelChanges")