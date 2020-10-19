args=commandArgs(trailingOnly=TRUE)

intercept=as.matrix(read.csv(paste(args[3],"_intercept.csv",sep=''),sep=',',header=TRUE))
coef=as.matrix(read.csv(paste(args[3],"_beta.csv",sep=''),sep=',',header=TRUE))
test=as.matrix(read.csv(args[1],sep=',',header=FALSE))



testdata<-test[,3:ncol(test)-1]

nclass<-length(intercept)
nlist<-seq(1,nclass)

calmat<-matrix(,nrow=nrow(testdata),ncol=nclass+1)


pred<-function(row,coefficients,intercept){
	tot<-(sum(coefficients*row)+intercept)
	return(exp(tot))}
for (eachind in 1:nrow(testdata)){
	l<-vector("list",nclass)
	
	for (eachclass in nlist){
	
	l[eachclass]<-pred(testdata[eachind,],coef[eachclass,],intercept[eachclass])
	}
	totprob<-Reduce("+",l)
	pl<-vector("list",nclass)
	for (eachagain in nlist){
		
	pl[eachagain]<-as.numeric(l[eachagain])/totprob}
			
	predclass<-(which.max(pl))

	calmat[eachind,1]<-predclass
	plr<-(as.numeric(unlist(pl)))

	
	calmat[eachind,2:(nclass+1)]<-plr}





rmat<-cbind(test[,1],calmat)
names<-c('index','class')
for (i in nlist){
	name<-paste('class',as.character(i),sep='_')
	names<-c(names,name)}

colnames(rmat)<-names
write.table(x=rmat,file=paste(args[2],sep=''),col.names=TRUE,row.names=FALSE,sep=',')






