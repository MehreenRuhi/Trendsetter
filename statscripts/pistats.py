import sys
import os
fname=sys.argv[1]
cwd=sys.argv[2]

f=open(cwd+'/'+fname+'.parse.tmp','r')
f3=open(cwd+'/'+fname+'_pistats.tmp','a')
from collections import Counter
import sys
import numpy as np
sites=[]
data=[]

def hetfunc(datasamp):
    	dsamp2=[]
	
	for tt in datasamp:
		dsamp2.append(map(float,tt))	
			
	
		
	heter=0
	for yy1 in dsamp2:
			
                cou=sum(yy1)
                p=cou/float(len(datasamp[0]))
                heter+=(2*p)*(1-p)
		heter2=heter/len(datasamp)
	
        f3.write(str(heter2))
        f3.write(',')
        f3.write('\n')
	
for i in f:
	sites.append(i.split(',')[0])
	data.append(i.rstrip('\n').split(',')[1])

ranges=range(0,len(sites),5)

for yy in range(0,(len(sites)/5)):
	
	d2=data[yy*5:(yy*5)+11]
	
	
	diversity=hetfunc(d2)


	
f.close()
f3.close()


