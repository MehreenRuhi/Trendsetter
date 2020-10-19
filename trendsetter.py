import os,sys
import argparse
import subprocess
parser = argparse.ArgumentParser(description='Trendsetter')
subparsers = parser.add_subparsers(help='sub-command help')
parser_calc = subparsers.add_parser('calcstat', help='calculate summary statistics from simulated data')
parser_calc.add_argument('fname', help='Simulated data file')
parser_calc.set_defaults(mode='calcstat')

parser_calc2 = subparsers.add_parser('calcstat_emp', help='calculate summary statistics from empirical data')
parser_calc2.add_argument('fname', help='data file')
parser_calc2.set_defaults(mode='calcstat_emp')

parser_train = subparsers.add_parser('train', help='Train classifier')
parser_train.add_argument('classifier_name',help='name of classifier')
parser_train.add_argument('trend_penalty',help='constant or linear')
parser_train.add_argument('num_stats',help='the number of different statistics (m)')
parser_train.add_argument('class_k',nargs='*',help='classes to be differentiated between (class_1.stats,class2.stats....class_k.stats)')
parser_train.set_defaults(mode='train')

parser_traincal=subparsers.add_parser('cal',help='platt scale classifier')
parser_traincal.add_argument('classifier_name',help='name of classifier')
parser_traincal.add_argument('class_k_cal',nargs='*',help='data for calibration (same order as training data)')
parser_traincal.set_defaults(mode='cal')

parser_test = subparsers.add_parser('test', help='test classifier')
parser_test.add_argument('test_file',help='summary statistics of testing dataset')
parser_test.add_argument('classifier_name',help='name of classifier')
parser_test.add_argument('testout',help='name of file to output results')
parser_test.set_defaults(mode='test')

parser_testcal = subparsers.add_parser('testcal', help='test classifier')
parser_testcal.add_argument('test_file',help='outputfile from test function (testout)')
parser_testcal.add_argument('classifier_name',help='name of classifier')
parser_testcal.set_defaults(mode='testcal')


ins = parser.parse_args()
allins = vars(ins)

if allins['mode'] == 'calcstat':
	
	calcinput = [allins['fname']]
	calculate = "python extrascripts/extend.py " + " ".join([str(c) for c in calcinput])
	subprocess.call(calculate.split())

elif allins['mode']=='calcstat_emp':
	cwd=os.getcwd()
	fname=allins['fname']
	calcin=[allins['fname'],cwd]
	parse='python statscripts/vcfparse.py '+ " ".join([str(c) for c in calcin])
	subprocess.call(parse.split())
	hstat='python statscripts/stats.py '+ " ".join([str(c) for c in calcin])
	subprocess.call(hstat.split())
	pistat='python statscripts/pistats.py '+ " ".join([str(c) for c in calcin])
	subprocess.call(pistat.split())
	rstat='python statscripts/rsqstats.py '+ " ".join([str(c) for c in calcin])
	subprocess.call(rstat.split())
	formstat='python statscripts/formtest.py '+ " ".join([str(c) for c in calcin])
	subprocess.call(formstat.split())
	os.remove(cwd+'/'+fname+'_rsqstats.tmp')
	os.remove(cwd+'/'+fname+'_pistats.tmp')
	os.remove(cwd+'/'+fname+'_hapstats.tmp')
	os.remove(cwd+'/'+fname+'.parse.tmp')
	os.remove(cwd+'/'+fname+'_sites.tmp')



elif allins['mode']=='train':
	fileslist= allins['class_k']
	numclass=str(len(fileslist))
	classifiername=allins['classifier_name']
	numstats=allins['num_stats']
	with open(classifiername+'_trainingdata','w') as trainingstats:
		for eachclassfile in fileslist:
			with open(eachclassfile) as currentfile:
				trainingstats.write(currentfile.read())
	standardize='python extrascripts/std.py '+classifiername +' '+numclass
	subprocess.call(standardize.split())
	inps=[classifiername,numstats]
	
	trainconstant=''	
	if allins['trend_penalty']=='constant':
		
		trainconstant="Rscript extrascripts/crossvalidategraph.R 1 "+" ".join([str(c) for c in inps])
		subprocess.call(trainconstant.split())
	elif allins['trend_penalty']=='linear':
		trainconstant="Rscript extrascripts/crossvalidategraph.R 2 "+" ".join([str(c) for c in inps])
		subprocess.call(trainconstant.split())
	os.remove(classifiername+'_trainingdata')
	os.remove(classifiername+'.data')
	
elif allins['mode']=='cal':
	fileslist= allins['class_k_cal']
	numclass=str(len(fileslist))
	cname=allins['classifier_name']
	calprobfile=open(cname+'cal.probs','a')	
	
	for eachtest in fileslist:
		tfile=eachtest.split('/')[-1]
		testout=tfile+'.probs'
		standardizetest=''
		standardizetest='python extrascripts/stdtest.py '+ cname+ ' '+eachtest
		subprocess.call(standardizetest.split())
		inps=[tfile,testout,cname]
		predict='Rscript extrascripts/predict.R std_' + " ".join([str(c) for c in inps])
		subprocess.call(predict.split())
		os.remove('std_'+tfile)
		with open(testout) as currentfile:
			
			currentitems=currentfile.readlines()[1:]
			for eachitem in currentitems:
				calprobfile.write("%s" % eachitem)
		os.remove(tfile+'.probs')		
	calprobfile.close()
	calclass='Rscript extrascripts/cal.R '+cname+' '+numclass
	subprocess.call(calclass.split())
	os.remove(cname+'cal.probs')
			
elif allins['mode']=='test':
	testfname=allins['test_file']
	testout=allins['testout']
	cname=allins['classifier_name']
	tfile=testfname.split('/')[-1]
	standardizetest=''
	inps=[tfile,testout,cname]
	predict=''
	standardizetest='python extrascripts/stdtest.py '+ cname+ ' '+testfname
	predict='Rscript extrascripts/predict.R std_' + " ".join([str(c) for c in inps])
	subprocess.call(standardizetest.split())
	subprocess.call(predict.split())
	os.remove('std_'+tfile)


elif allins['mode']=='testcal':
	testfname=allins['test_file']
	cname=allins['classifier_name']
	tfile=testfname.split('/')[-1]
	predict=''
	predict='Rscript extrascripts/predcal.R '+testfname+' '+cname
	subprocess.call(predict.split())
	


