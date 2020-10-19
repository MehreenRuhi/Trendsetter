from collections import Counter
import sys


def hetfunc(datasamp):
	heter=0
        for tt in range(len(datasamp[0])):
 	       cou=0
               for ttt in range(len(datasamp)):
			cou+=float(datasamp[ttt][tt])
               p=cou/float(len(datasamp))


               heter+=(2*p)*(1-p)
	return heter

def rsquared(midsnp,datasamp):

        all1=0.00
        all2=0.00
        allb1=0.00
        allb2=0.00
        a1b1=0.00
                 
        for ll in range(len(datasamp)):
               	if datasamp[ll]=='0':
		        all1+=1.000
                if midsnp[ll]=='0':
                        allb1+=1.000
                if datasamp[ll]=='0' and midsnp[ll]=='0':
                        a1b1+=1.000
	allb2=float(len(datasamp))-allb1
        all2=float(len(datasamp))-all1
        fh= a1b1/len(datasamp)
        sh=(all1/len(datasamp))*(allb1/len(datasamp))
        D=fh-sh
        Dsq=D*D
        bot=((all1/len(datasamp))*(allb1/len(datasamp))*(all2/len(datasamp))*(allb2/len(datasamp)))
        rsq=Dsq/bot
	return rsq                     

file1=sys.argv[1]

f2=open(file1+'.stats','a')


f=open(file1,'r')

g=f.readlines()
pop=int(g[0].split(' ')[1])+2
for i in range(len(g)):
	segs1=[]
	data=[]
	if g[i][0:3]=='seg':
		
		segs1=g[i+1].rstrip(' \n').lstrip('positions: ').split(' ')
		segs=map(float,segs1)
		seglistlen=len(segs)
		if seglistlen>1010:
			mid= seglistlen/2
			data=g[i+2:i+pop]
			
			ranges=range(-505,510,5)
			f2.write(str(segs1[mid]))
			f2.write(',')
			for win in range(201):
				midsnp=[]
				winsegs=[]
				for eachline in data:
					midsnp.append(eachline[mid])
					winsegs.append(eachline[mid+ranges[win]:(mid+ranges[win+2])])
				rsqlist=[]
				for eachwinseg in range(len(winsegs[0])):
					eachcolwin=[]
					for eachwin in winsegs:
					
						eachcolwin.append(eachwin[eachwinseg])
					rsqlist.append(rsquared(midsnp,eachcolwin))
				pi=hetfunc(winsegs)
				avepi=pi/len(winsegs[0])
				f2.write(str(avepi))
				f2.write(',')
				rsq= sum(rsqlist)/len(winsegs[0])	
				f2.write(str(rsq))
				f2.write(',')	
		                count=Counter(winsegs)
                                counts=[]

                                for l in set(winsegs):
                                        counts.append(count[l])
                                f2.write(str(len(set(winsegs))))
                                f2.write(',')
                                sortedcount= sorted(counts,reverse=True)
                                hlist=[]
                                for each in sortedcount:
                                        hlist.append(each)
				hlist2=[eachH/(float(len(winsegs))) for eachH in hlist]
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
		else:
			pass
		

f.close()	
f2.close()
