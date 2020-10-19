args=commandArgs(trailingOnly=TRUE)
library(glmnet)
testdata=as.data.frame(read.csv(args[1],sep=',',header=TRUE))

cal<-readRDS(paste(args[2],'calmodel.rds',sep=''))
newdata=as.matrix(testdata[,3:5])

newdataprobs<-as.data.frame(predict(cal,newdata,type='response'))

newprobs<-cbind(testdata[,1],newdataprobs)
colnames(newprobs)<-c('index','class1','class2','class3')
write.table(newprobs,file=paste(args[1],'.calprobs',sep=''),sep=',',col.names=TRUE,row.names=FALSE)
