workingdir=getwd()
if("fusedlasso" %in% rownames(installed.packages()) == FALSE) {install.packages(paste(workingdir,'/fusedlasso-master',sep=''),dependencies=TRUE,repo=NULL)}

args = commandArgs(trailingOnly = TRUE)
infile <- args[2]
d <- as.integer(args[1])
difstats<-as.integer(args[3])


foldCV = 10  # To make 10-fold cross validation

if(d != 1 & d != 2) {
    print("ERROR: the trend-filtering derivative should be of order d = 1 or d = 2")
    q()
}



predict<-function(row,coefficients,intercept){
	tot<-(sum(coefficients*row)+intercept)
	return(exp(tot))}
	

x<-as.matrix(read.csv(paste(infile,'.data',sep=''),sep=',',header=FALSE))

library(fusedlasso)


set.seed(1)

# Assign IDs for cross validation
classSet <- unique(x[,1]) # Get the set of unique class labels
id <- c()
for(i in classSet) {
    myPerm <- sample(which(x[,1] == i))    
    myChoice <- rep(1:foldCV,length(myPerm))  
    
    sizeCV <- floor(length(myPerm)/foldCV)
    myChoice <- rep(foldCV,length(myPerm))
    for(j in 1:foldCV) {
        for(k in 1:sizeCV) {
            myChoice[(j-1)*sizeCV + k] = j
        }
    }

    for(j in 1:length(myPerm)) {
        id[myPerm[j]] = myChoice[j] 
    }
}

lambdalist<-c(0.0001,0.0005, 0.001,0.005,0.01,0.05)
lambdalist2<-c( 0.0001,0.0005,0.001,0.005,0.01,0.05)
lambdalistExpanded <- c()
lambdalist2Expanded <- c()
numLambda <- 0
for (i in 1:length(lambdalist)){
	for (j in i:length(lambdalist2)){
        numLambda<<-numLambda+1
        lambdalistExpanded[numLambda] <- lambdalist[i]
        lambdalist2Expanded[numLambda] <- lambdalist2[j]
    }
}

resmat<-matrix(,nrow=numLambda,ncol=3)
resmat[is.na(resmat)]<-0

groupposition<-c(rep(1:difstats,(ncol(x)-1)/difstats))

library(foreach)
library(doParallel)

numCores <- detectCores() - 1   # GET MAXIMUM NUMBER OF AVAILABLE CORES
myClust <- makeCluster(numCores)    # INITIATE A CLUSTER
registerDoParallel(myClust)

traindata<-x[,2:ncol(x)]
trainclass<-x[,1]
# Run the CV across lambda value pairs in parallel, and then store the results by rows into the resultant matrix
resmat <- foreach(i=seq(1:numLambda), 
                  .combine = rbind,
                  .packages = c("fusedlasso")) %dopar%
{
    marginalResult <- matrix(0,nrow=1,ncol=3)
        
	for (cross in 1:foldCV) {
		num<-which (!(id %in% cross))
		notnum<-which(id %in% cross)
		traindata<-x[num,2:ncol(x)]
		trainclass<-x[num,1]
		testdata<-x[notnum,2:ncol(x)]
		testclass<-x[notnum,1]
	    
        res<-fusedlasso(traindata,trainclass,lambda.lasso=lambdalistExpanded[i],lambda.fused=lambdalist2Expanded[i],groups=groupposition,family="multinomial",d=d,gamma1=1,gamma2=1)
        		
xcount<-0

for (eachind in 1:nrow(testdata)){
	totprob<-0
	l<-vector("list",length(classSet))
	
	for (eachclass in classSet){
		l[eachclass+1]<-predict(testdata[eachind,],res$beta[eachclass+1,],res$intercept[eachclass+1])
		
			
		}
			totprob<-Reduce("+",l)
			
			pl<-vector("list",length(classSet))
	for (eachagain in classSet){
		
		pl[eachagain+1]<-as.numeric(l[eachagain+1])/totprob}
				
		predclass<-(which.max(pl)-1)
		actualclass<-(testclass[eachind])
		
		if (predclass==actualclass){
			xcount <- xcount+1
            }
    	}	
        
        marginalResult[1] <- lambdalistExpanded[i]
		marginalResult[2] <- lambdalist2Expanded[i]
		marginalResult[3] <- marginalResult[3] + xcount
    }

    marginalResult   # Return result for this lambda pair to be stored in combined result matrix
}

stopCluster(myClust)  # STOP THE CLUSTER

bestind<-(which.max(resmat[,3]))

print(paste("Identified best pair of lambda values (", resmat[bestind,1], ",", resmat[bestind,2], ")",sep=""))
print("Computing model based on these lambda values")

traindata<-x[,2:ncol(x)]
trainclass<-x[,1]
bestres<-fusedlasso(traindata,trainclass,lambda.lasso=resmat[bestind,1],lambda.fused=resmat[bestind,2],groups=groupposition,family="multinomial",d=d,gamma1=1,gamma2=1)

write.table(bestres$beta, file=paste(infile, "_beta.csv",sep=""),sep=',',row.names=T,col.names=T)
write.table(bestres$intercept,file=paste(infile, "_intercept.csv",sep=""),sep=',',row.names=T,col.names=T)



