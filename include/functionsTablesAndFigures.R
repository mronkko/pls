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
writeComparisonTable <- function(data,variables,file,analysisTypes,labels){

	tableData<-NULL

	tempdata<-data[,c("designNumber","analysis",variables)]
	tempdata<-aggregate(tempdata, by=list(data$designNumber,data$analysis),  FUN=mean, na.rm=TRUE)
	
	comparisonData<-reshape(tempdata[,c("designNumber","analysis",variables)],v.names=variables,idvar="designNumber",timevar="analysis",direction="wide")
	
		
	for(i in 1:length(analysisTypes)){
		tableRow=data.frame(rownames=labels[[analysisTypes[i]]])
		for(j in 1:length(variables)){
			
			varname<-variables[j]
			
			print(paste(varname,i,sep="."))
			print(summary(comparisonData))
			print(comparisonData[,paste(varname,i,sep=".")])
			tableRow<-cbind(tableRow,t(quantile(comparisonData[,paste(varname,i,sep=".")],probs=c(0.05,0.5,0.95),na.rm=TRUE)))
		

			# Calculate how often this analysis results in larger results than every other analysis

			tableRow<-cbind(tableRow,sum(comparisonData[,paste(varname,i,sep=".")] == apply(comparisonData[,paste(varname,1:length(analysisTypes),sep=".")], 1, max)))

			# Calculate how often this analysis results in smaller results than every other analysis

			tableRow<-cbind(tableRow,sum(comparisonData[,paste(varname,i,sep=".")] == apply(comparisonData[,paste(varname,1:length(analysisTypes),sep=".")], 1, min)))

	
		}
		tableData<-rbind(tableData,tableRow)
	}
	print(xtable(tableData),file=paste("results/",file,"_full.tex",sep=""))
	print(xtable(tableData),file=paste("results/",file,"_body.tex",sep=""),include.rownames=FALSE,
hline.after=NULL,only.contents=TRUE,include.colnames=FALSE)
	return(tableData)

}