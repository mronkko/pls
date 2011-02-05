#
# Functions used when to drawing tables and figures
#

library(xtable)

writeComparisonTable <- function(data,variables,file,analysisTypes){

	tableData<-NULL

	for(i in 1:length(analysisTypes)){
		tableRow=data.frame(rownames=analysisTypes[i])
		thisData<-fullData[fullData$analysis==i,]
		for(j in 1:length(variables)){
			varname<-variables[j]
			tableRow<-cbind(tableRow,t(quantile(thisData[,varname],probs=c(0.05,0.5,0.95))))
			for(k in 1:length(analysisTypes)){
				print(paste("Calculating table. Indexes are",i,j,k))
				if(k!=i){
					
					# Calculate how often this analysis results in larger results than every other analysis
					
					comparisonData<-fullData[fullData$analysis==k,]
					mergedData<-merge(thisData,comparisonData,by =c("replication","designNumber","construct"))
					comparison<-mergedData[,paste(varname,".x",sep="")]>mergedData[,paste(varname,".y",sep="")]
					tableRow<-cbind(tableRow,sum(comparison,na.rm=TRUE)/sum(!is.na(comparison)))
				}
				else{
					tableRow<-cbind(tableRow,NA)
				}
				colnames(tableRow)[ncol(tableRow)]<-analysisTypes[k]
	
			}
		}
		tableData<-rbind(tableData,tableRow)
	}
	print(xtable(tableData),file=paste("results/",file,"_full.tex",sep=""))
	print(xtable(tableData),file=paste("results/",file,"_body.tex",sep=""),only.contents=TRUE)
	return(tableData)

}