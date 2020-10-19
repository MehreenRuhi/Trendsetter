
import sys
import os
fname=sys.argv[1]
cwd=sys.argv[2]

f=open(cwd+'/'+fname+'_rsqstats.tmp','r')
g=open(cwd+'/'+fname+'_pistats.tmp','r')
h=open(cwd+'/'+fname+'_hapstats.tmp','r')
site=open(cwd+'/'+fname+'_sites.tmp','r')
w=open(fname+'.stats','a')

pilist=[]
sitelist=[]
haplist=[]
rsqlist=[]
for i in f:
	rsqlist.append(i.rstrip(',\n').split(','))

for j in g:
	pilist.append(j.rstrip('\n'))


for k in h:
	haplist.append(k.rstrip('\n'))
for r in site:
	sitelist.append(r.rstrip('\n'))

for t in range(len(rsqlist)-101):	
	pileft=[]
	rsqleft=[]
	haps=[]
	rsqlistwhole=[]
	haps.append(haplist[t:t+201])
	piwhole=[]
	piwhole.append(pilist[t:t+201])
	        	
	for it2 in range(1,101):
		
                rsqleft.append(rsqlist[t-1+it2][-it2])
		
	rsqlistwhole=rsqleft+rsqlist[t+101]
		
	w.write(sitelist[t+101])
	w.write(',')
	for d in range(len(haps[0])):
		w.write(piwhole[0][d])
		w.write(rsqlistwhole[d])
		w.write(',')
		w.write(haps[0][d])
	w.write('\n')	


site.close()
f.close()
w.close()
g.close()
h.close()
