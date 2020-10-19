import sys

import numpy
fname=sys.argv[1]
nclass=int(sys.argv[2])
f=open(fname+'_trainingdata','r')
g=open(fname+'.std','a')
overall=[]
for i in f:
	
	if i[0]!='n':
		x=i.rstrip(',\n')
		y=map(float,x.split(',')[1:])
		
		overall.append(y)

s=numpy.array(overall)
t=numpy.std(s,axis=0)
u=numpy.mean(s,axis=0)		
v=s-u
w=v/t


clist=range(0,nclass)
classes=[]
for cl in clist:
	classes+=([cl]*(len(overall)/nclass))

wx=numpy.column_stack([classes,w])

numpy.savetxt(g,t)
numpy.savetxt(g,u)
numpy.savetxt(fname+'.data',wx, delimiter=',')
