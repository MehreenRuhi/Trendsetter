import sys
import os
cwd=sys.argv[2]

fname=sys.argv[1]
out=open(cwd+'/'+fname+'.parse.tmp','a')
with open(fname,'r') as i:
	for j in i:
		popinfor=[]
		if j[0]!='#':
			geneline=j.rstrip('\n').split('\t')[9:]
			for eachind in geneline:
				for hap in eachind.split('|'):
					popinfor.append(hap)
		if len(set(popinfor))==2:
			out.write(str(j.split('\t')[1]))
			out.write(',')
			out.write(''.join(popinfor))
			out.write('\n')
out.close()
			
