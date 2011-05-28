#
# Functions used when to drawing tables and figures
#

library(xtable)

writeDescriptivesTable <- function(data,variables,file,analysisTypes,labels){

	tableData<-NULL
	probs=c(.05 , .5, .95)
	
	command<-NULL
	pos<-NULL
	
	# Each variable gets a block of the table. The table in total has two "mega columns"
	
	tableTwoBlocks<-NULL
	
	for(i in 1:length(variables)){
		tableRow<-NULL
		varname<-variables[i]
		for(j in 1:length(analysisTypes)){
			thisData<-data[data$analysis==j,varname]
			tableRow<-cbind(tableRow,quantile(thisData,probs=probs))
		}
		colnames(tableRow)<-analysisTypes
		tableRow<-data.frame(paste(probs*100,"%",sep=""),tableRow)

		if(i %% 2== 0){
			tableData<-rbind(tableData,cbind(tableOddBlock,tableRow))
			pos=c(pos,(i/2-1)*3)
			command=c(command,paste("\\multicolumn{",ncol(tableRow),"}{l}{",labels[[variables[i-1]]],"}&\\multicolumn{",ncol(tableRow),"}{@{}l}{",labels[[varname]],"}\\\\"))
		}
		else{
			tableOddBlock<-tableRow
		}
	}
	
	#If there are odd number of things to include, there is one "left over part"

	if(length(variables) %% 2 ==1){
		nullBlock<-matrix(rep(NA,ncol(tableRow)*nrow(tableRow)),ncol=ncol(tableRow))
		
		print(nullBlock)
		
		colnames(nullBlock)<-colnames(tableRow)
		
		tableData<-rbind(tableData,cbind(tableRow,nullBlock))
		pos=c(pos,((i+1)/2-1)*3)
		command=c(command,paste("\\multicolumn{",ncol(tableRow),"}{l}{",labels[[varname]],"}\\\\"))
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
			
#			print(paste(varname,i,sep="."))
#			print(summary(comparisonData))
#			print(comparisonData[,paste(varname,i,sep=".")])
#			print(tableRow)
			print(t(quantile(comparisonData[,paste(varname,i,sep=".")],probs=c(0.05,0.5,0.95),na.rm=TRUE)))
			
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

#
# Strips lmer results to bare minimum that is required for drawing regression tables.
#

setClass("strippedlmer",representation=representation(ll="logLik",coefs="matrix",REmat="matrix",resid="numeric"))

striplmer <- function(x){
	return(new("strippedlmer",ll=logLik(x),coefs=summary(x)@coefs,REmat=summary(x)@REmat,resid=x@resid))
}
#
# Two helper functions relating to stripping lmer objects
#

summary.strippedlmer <- function(x){
	return(x)
}

logLik.strippedlmer <- function(x){
	return(x@ll)
}

#
# Combines results from lmer 
#
# Source:
# http://www.rensenieuwenhuis.nl/r-sessions-31-combining-lmer-output-in-a-single-table/
#

combine.output.lmer <- function(models, labels=FALSE)
{

fix.coef <- lapply(models, function(x) summary(x)@coefs)
var.coef <- lapply(models, function(x) summary(x)@REmat)
n.par <- dim(summary(models[[1]])@coefs)[2]

ifelse(labels==FALSE,
fix.labels <- colnames(summary(models[[1]])@coefs),
fix.labels <- labels)

var.labels <- colnames(var.coef[[1]])

# Creating table with fixed parameters
output.coefs <- data.frame(Row.names=row.names(fix.coef[[1]]))
for (i in 1:length(models))
{

a <- fix.coef[[i]]
colnames(a) <- paste("Model", i, fix.labels)
output.coefs <- merge(output.coefs, a, by.x=1, by.y=0, all=T, sort=FALSE)

}
output.coefs[,1] <- as.character(output.coefs[,1])
output.coefs[dim(output.coefs)[1]+2, 1] <- "Loglikelihood"
LL <- unlist(lapply(models, function(x) as.numeric(logLik(x))))
output.coefs[dim(output.coefs)[1], 1:length(models)*n.par-n.par+2] <- LL

# Creating table with random parameters
output.vars <- data.frame(var.coef[[1]])[,1:2]
for (i in 1:length(models))
{

a <- var.coef[[i]]
colnames(a) <- paste("Model", i, var.labels)
output.vars <- merge(output.vars, a, by.x=1:2, by.y=1:2, all=T, sort=FALSE)

}

# Combining output.coefs and output.vars
n.cols <- dim(output.coefs)[2]
n.coefs <- dim(output.coefs)[1]
n.vars <- dim(output.vars)[1]

output <- matrix(ncol=n.cols +1 , nrow=n.vars+n.coefs+2)

output[1:n.coefs, -2] <- as.matrix(output.coefs)
output[n.coefs+2, 1] <- "Variance Components"
output[(n.coefs+3) : (n.coefs+n.vars+2), 1:2] <- as.matrix(output.vars[,1:2])
output[
(n.coefs+3) : (n.coefs+n.vars+2),
which(rep(c(1,1,rep(0, n.par-2)),length(models))!=0)+2] <- as.matrix(output.vars[,c(-1,-2)])

colnames(output) <- c("Parameter", "Random", colnames(output.coefs)[-1])

return(output)
}
