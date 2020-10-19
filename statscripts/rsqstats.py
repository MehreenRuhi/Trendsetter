import sys
import os
fname=sys.argv[1]
cwd=sys.argv[2]

f=open(cwd+'/'+fname+'.parse.tmp','r')
f3=open(cwd+'/'+fname+'_rsqstats.tmp','a')

from collections import Counter
import sys
import numpy as np
sites=[]
data=[]
	

def rsquared2(datasamp):
        all1=0.00
        all2=0.00
        allb1=0.00
        allb2=0.00
        a1b1=0.00
	d3h=''
        for ll in datasamp:
		d3h+=ll[-1]
                if ll[0]=='0':
                        all1+=1.000
                if ll[-1]=='0':
                        allb1+=1.000
                if ll[0]=='0' and ll[-1]=='0':
                        a1b1+=1.000
        print all1, allb1, a1b1
	
	allb2=float(len(datasamp))-allb1
        all2=float(len(datasamp))-all1
        fh= a1b1/len(datasamp)
        sh=(all1/len(datasamp))*(allb1/len(datasamp))
        D=fh-sh
        Dsq=D*D
        bot=((all1/len(datasamp))*(allb1/len(datasamp))*(all2/len(datasamp))*(allb2/len(datasamp)))
        rsq=Dsq/bot
        print len(datasamp),'datasamp'
	return rsq
	


def rsquared(d2,d3):
	all1=0.00
        all2=0.00
        allb1=0.00
        allb2=0.00
        a1b1=0.00
	for ll in range(len(d2)):
		if d2[ll]=='0':
			all1+=1.000
		if d3[ll]=='0':
			allb1+=1.000
		if d2[ll]=='0' and d3[ll]=='0':
			a1b1+=1.000
 
        allb2=float(len(d2))-allb1
        all2=float(len(d2))-all1
        fh= a1b1/len(d2)
        sh=(all1/len(d2))*(allb1/len(d2))
        D=fh-sh
        Dsq=D*D
        bot=((all1/len(d2))*(allb1/len(d2))*(all2/len(d2))*(allb2/len(d2)))
        rsq=Dsq/bot
	
        return rsq


for i in f:
	sites.append(i.split(',')[0])
	data.append(i.rstrip('\n').split(',')[1])

ranges3=range(0,101)

#for yy in range(0,(len(sites)/5)-500):
for yy in range(0,(len(sites)-500)/5):
	d1=data[yy*5]

	for yyy in ranges3:
		
		groupdata=data[(yy*5)+(yyy*5): ((yy*5)+(yyy*5)+11)]
		rsqlist=[]
		for eachsnp in groupdata:
			d3=eachsnp
			indrsq=rsquared(d1,d3)
			rsqlist.append(indrsq)
		rsqave=sum(rsqlist)/len(groupdata)
		f3.write(str(rsqave))
		f3.write(',')
	f3.write('\n')
	
f.close()
f3.close()


