
import sys

import os
fname=sys.argv[1]
cwd=sys.argv[2]

f=open(cwd+'/'+fname+'.parse.tmp','r')
f2=open(cwd+'/'+fname+'_hapstats.tmp','a')
f3=open(cwd+'/'+fname+'_sites.tmp','a')
from collections import Counter

import numpy as np
sites=[]
data=[]

for i in f:
	sites.append(i.split(',')[0])
	data.append(i.rstrip('\n').split(',')[1])


for t in range(0,(len(sites)/5)):

	f3.write(sites[t*5])
	f3.write('\n')	
	
	d= data[t*5:(t*5)+11]

	d=[map(None, dd) for dd in d]

	hapmat=np.matrix(d)
	hapmat2=hapmat.T
	haplis=np.array(hapmat2).tolist()
	haplist=[''.join(y) for y in haplis]
	
		
	count=Counter(haplist)
        counts=[]

        for l in set(haplist):
        	counts.append(count[l])
        f2.write(str(len(set(haplist))))
        f2.write(',')
        sortedcount= sorted(counts,reverse=True)
        hlist=[]
        for each in sortedcount:
        	hlist.append(each)
        hlist2=[eachH/(float(len(haplist))) for eachH in hlist]
        hsqlist=[j ** 2 for j in hlist2]

        h2sq=hsqlist[1:]
        h2sum=float(sum(h2sq))
        h1=float(sum(hsqlist))
        h21=h2sum/h1
        one=hlist2[0]
        two=hlist2[1]
        h12sqlist=h2sq[1:]
        h12part1=(one+two)*(one+two)
        h12part2=sum(h12sqlist)
        h12=h12part1+h12part2
        f2.write(str(h1))
        f2.write(',')
        f2.write(str(h12))
        f2.write(',')
        f2.write(str(h21))
        f2.write(',')
	f2.write('\n')

	
f.close()
f2.close()
f3.close()


