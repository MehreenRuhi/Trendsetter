args=commandArgs(trailingOnly=TRUE)
library(glmnet)
caldata=as.data.frame(read.csv(paste(args[1],"cal.probs",sep=''),sep=',',header=FALSE))
numclass=args[2]
nlist<-seq(1,numclass)
rowscal<-nrow(caldata)
y<-rep(nlist, each=(rowscal/as.numeric(numclass)))
x<-as.matrix(caldata[,3:5])

cal=glmnet(x,y,family="multinomial",lambda=0)

saveRDS(cal,paste(args[1],'calmodel.rds',sep=''))
