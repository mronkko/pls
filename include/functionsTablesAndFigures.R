#
# Functions used when to drawing tables and figures
#

library(xtable)

writeDescriptivesTable <- function(data,variables,file,analysisTypes,labels){

	tableData<-NULL
	probs=c(.05 , .5, .95)
	
	command<-NULL
	pos<-NULL
	
	# Each variable gets a block of the table
	
	for(i in 1:length(variables)){
		tableRow<-NULL
		varname<-variables[i]
		for(j in 1:length(analysisTypes)){
			thisData<-data[data$analysis==j,varname]
			tableRow<-cbind(tableRow,quantile(thisData,probs=probs))
		}
		colnames(tableRow)<-analysisTypes
		tableRow<-data.frame(paste(probs*100,"%",sep=""),tableRow)
		pos=c(pos,(i-1)*3)
		command=c(command,paste("\\multicolumn{",ncol(tableRow),"}{l}{",labels[[varname]],"}\\\\"))
		tableData<-rbind(tableData,tableRow)
	}
	
	add.to.row<-list(pos=as.list(pos),command=command)
	
	print(xtable(tableData,digits=3),file=paste("results/",file,"_full.tex",sep=""),add.to.row=add.to.row)
	print(xtable(tableData,digits=3),file=paste("results/",file,"_body.tex",sep=""),include.rownames=FALSE,
hline.after=NULL,only.contents=TRUE,include.colnames=FALSE,add.to.row=add.to.row)
	return(tableData)

}
writeComparisonTable <- function(data,variables,file,analysisTypes){

	doOnDesignLevel<-TRUE
	tableData<-NULL

	for(i in 1:length(analysisTypes)){
		tableRow=data.frame(rownames=analysisTypes[i])
		thisData<-data[data$analysis==i,]
		for(j in 1:length(variables)){
			
			varname<-variables[j]
			tableRow<-cbind(tableRow,t(quantile(thisData[,varname],probs=c(0.05,0.5,0.95),na.rm=TRUE)))
			for(k in 1:length(analysisTypes)){
				print(paste("Calculating table. Indexes are",i,j,k))
				if(k!=i){
					
					# Calculate how often this analysis results in larger results than every other analysis

					comparisonData<-data[data$analysis==k,]
					mergedData<-merge(thisData,comparisonData,by=intersect(names(thisData),c("replication","designNumber","construct","to","from")))
					comparison<-mergedData[,paste(varname,".x",sep="")]>mergedData[,paste(varname,".y",sep="")]
					
					if(doOnDesignLevel){
						# in how many designs this is better than the alternative
						agg<-aggregate(comparison, by=list(mergedData$designNumber),  FUN=mean, na.rm=TRUE)
						tableRow<-cbind(tableRow,sum(agg[,2]>.5))
					}
					else{
						tableRow<-cbind(tableRow,sum(comparison,na.rm=TRUE)/sum(!is.na(comparison)))
					}
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
	print(xtable(tableData),file=paste("results/",file,"_body.tex",sep=""),include.rownames=FALSE,
hline.after=NULL,only.contents=TRUE,include.colnames=FALSE)
	return(tableData)

}