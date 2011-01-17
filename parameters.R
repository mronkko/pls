# Overall study related parameters

replications<-1

# Model related parameters ( 3 x 3 x 3 = 27 combinations)

numberOfConstructs <-c(4,8,12)
expectedNumberOfOutgoingPaths <-c(1,2,3)
populationPathValues<-list(c(-.5,.5),c(0,.5),c(0,0))


# Data related parameters (4 combinations)

sampleSizes<-c(25,50,100,250)

# Measurement model related parameters. These all are population values. ( 3 x 3 x 3 x 3 x 3 = 243 combinations)

indicatorCounts<-c(3,4,6)
factorLoadings<-c(.4,.6,.8)
factorLoadingIntervals<-c(0,.1,.2)
maxErrorCorrelations<-c(0,.2,.4)
methodVariances<-c(0,.1,.2)

# Parameters related to tested models (3 x 3 = 9 combinations)

omittedPathsShare<-c(0,.25,.5)
extraPaths<-c(1,2,3)
