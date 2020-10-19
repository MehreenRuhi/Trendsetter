
import sys
import numpy

tin=sys.argv[2]
fname=sys.argv[1]
outfname=tin.split('/')[-1]
f=open(tin,'r')
g=open(fname+'.std','r')
sites=[]
overall=[]

	

for i in f:
	x=i.strip(',\n')
	sites.append(x.split(',')[0])
	y=map(float,x.split(',')[1:])
	overall.append(y)

stdlist=[]
for j in g:
	jk=j.rstrip('\n')
	jy=float(jk)
	stdlist.append(jy)
	

t=numpy.array(stdlist[:(len(stdlist)/2)])
u=numpy.array(stdlist[(len(stdlist)/2):])
s=numpy.array(overall)

v=s-u
w=v/t

sitearray= numpy.matrix(sites)
with open('std_'+outfname,'a') as writefile:
	for eachrow in range(len(sites)):
		writefile.write(sites[eachrow])
		writefile.write(',')
		for eachentry in range(len(w[0])):
			
			writefile.write(str(w[eachrow][eachentry]))
			writefile.write(',')
		writefile.write('\n')
	

	

f.close()
g.close()

