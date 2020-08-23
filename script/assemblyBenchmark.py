__author__		= "Chayan Kumar Saha"
__copyright__	= "X-genomics"
__email__		= "chayan.sust7@gmail.com"

'''
Tools Required to be installed:

SraToolkit
FastQC
multiQC
Trimmomatic



'''


import argparse
import glob
import os, sys, os.path
import queue as Queue
import subprocess
from queue import Queue,Empty
import threading



usage= ''' Description: Download fastq files from NCBI using SRA number and proceeds with assembly benchmarking pipeline '''


parser = argparse.ArgumentParser(description=usage)
parser.add_argument("-i", "--input", required= True, help=" List of SRA files ")
parser.add_argument("-ap", "--adapterPE", required= True, help=" ILLUMINACLIP PE Adapter File ")
parser.add_argument("-as", "--adapterSE", required= True, help=" ILLUMINACLIP SE Adapter File ")
parser.add_argument("-c", "--cpu", help="Maximum number of parallel CPU workers to use for multithreads. ")
parser.add_argument("-k", "--keep", action="store_true", help=" If you want to keep the intermediate files use [-k]. By default it will remove. ")
parser.add_argument("-v", "--version", action="version", version='%(prog)s 1.0.0')
parser.add_argument("-vb", "--verbose", action="store_true", help=" Use this option to see the work progress for each query as stdout. ")
args = parser.parse_args()
parser.parse_args()

core=int(args.cpu)

def worker_func():
	while not stopped.is_set():
		try:
			# use the get_nowait() method for retrieving a queued item to
			# prevent the thread from blocking when the queue is empty
			com = q.get_nowait()
		except Empty:
			continue
		try:
			os.system(com)
			#subprocess.Popen(com, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
		except Exception as e:
			print('Error running command', str(e))
		finally:
			q.task_done()

'''
#Input File should be like this
SRR10971381	PAIRED
SRR11826840	PAIRED
SRR11851931	PAIRED
SRR11851968	PAIRED
SRR11859133	PAIRED
SRR12100425	SINGLE
SRR12100424	SINGLE
SRR12100423	SINGLE
'''

#sraSet=set() #Make a set of SRA IDs, duplicates are removed
sraSet={}
with open(args.input, 'r') as sraIn:
	for line in sraIn:
		if line[0]!='':
		#Line=line.rstrip()
			Line=line.rstrip().split('\t')
			#sraSet.add(Line)
			sraSet[Line[0]]=Line[1]


if args.verbose:
	print('\n>> '+args.input+ ' file contains '+str(len(sraSet))+ ' IDs \n')

#Creates a directory and Downloads fastq.gz files inside the directory, so each sample will be kept seprate in a folder

dloadSet=set()
makedataDir= './Data'
if not os.path.isdir(makedataDir):
	os.mkdir(makedataDir)
else:
	for items in (glob.glob(makedataDir+'/'+'*.gz')):
		if items:
			if items.split('/')[-1][-11]=='_':
				fileName=items.split('/')[-1].split('_')[0]
				dloadSet.add(fileName)
			else:
				fileName=items.split('/')[-1].replace('.fastq.gz','')
				dloadSet.add(fileName)
		else:
			fileName='Null'
			dloadSet.add(fileName)

scount=0
for sraIDs in sraSet:
	if sraIDs not in dloadSet:
		scount+=1
		dCommand="fastq-dump --split-3 %s --outdir %s --gzip" %(sraIDs,makedataDir)
		if args.verbose:
			print('\tDownloading ... '+sraIDs+' ('+str(scount)+'/'+str(len(sraSet))+')')
		dcom=subprocess.Popen(dCommand, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
		dcom.wait()
	else:
		if args.verbose:
			print('\t'+sraIDs+' downloaded previously')

print('\n>> Total '+str(scount)+' files were Downloaded from NCBI !! \n')

#sys.exit()
# FastQC run

prevFastQC=set()
makeFastQCDir= './Data_FastQC_report'
if not os.path.isdir(makeFastQCDir):
	os.mkdir(makeFastQCDir)
else:
	for items in (glob.glob(makeFastQCDir+'/'+'*html')):
		if items:
			fileName=items.split('/')[-1].split('_')[0]
			prevFastQC.add(fileName)
		else:
			fileName='Null'
			prevFastQC.add(fileName)

fastQC_com=[]


for sraIDs in sraSet:
	if sraIDs in prevFastQC:
		pass
	else:
		fastQC_Command="fastqc %s -o %s" %(makedataDir+'/'+sraIDs+'*', makeFastQCDir)
		fastQC_com.append(fastQC_Command)

print('\n>> FastQC Running for Data .. ')

if fastQC_com:
	thread_count = core # maximum parallel threads
	stopped = threading.Event()

	q = Queue()
	print('-- Processing QC with FastQC : '+ (str(len(fastQC_com))+ ' tasks in thread queue with '+ str(thread_count))+ ' thread limit')
	for item in fastQC_com:
		q.put(item)
	for x in range(thread_count):
		t = threading.Thread(target=worker_func)
		# t.daemon = True #Enable to run threads as daemons
		t.start()
	q.join()       # block until all tasks are done
	stopped.set()

# multiQC run

makeMultiQCInitialDir= './Data_multiQC_report'
if not os.path.isdir(makeMultiQCInitialDir):
	os.mkdir(makeMultiQCInitialDir)

multiQC_command_init="multiqc -f %s -o %s" %(makeFastQCDir, makeMultiQCInitialDir)
print('\n>> MultiQC Running for Data .. ')
multiQC_command_run=subprocess.Popen(multiQC_command_init, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
multiQC_command_run.wait()


#Trimmomatic
prevTrim=set()
makeTrimmedFastQDir= './TrimmedData'
if not os.path.isdir(makeTrimmedFastQDir):
	os.mkdir(makeTrimmedFastQDir)
else:
	for items in (glob.glob(makeTrimmedFastQDir+'/'+'*.gz')):
		if items:
			fileName=items.split('/')[-1].split('_')[0]
			prevTrim.add(fileName)
		else:
			fileName='Null'
			prevTrim.add(fileName)

fastQ_1="*_1.fastq.gz"
fastQ_2="*_2.fastq.gz"
fastQ_S=".fastq.gz"

print('\n>> Trimming Data with trimmomatic .. ')



for sraIDs in sraSet:
	if sraSet[sraIDs][0]=='P': #PAIRED
		if sraIDs in prevTrim:
			pass
		else:
			fastQ1=(glob.glob(makedataDir+'/'+sraIDs+fastQ_1))[0]
			fastQ2=(glob.glob(makedataDir+'/'+sraIDs+fastQ_2))[0]
			fpgz=makeTrimmedFastQDir+'/'+sraIDs+'_1P.fastq.gz'
			fugz=makeTrimmedFastQDir+'/'+sraIDs+'_1U.fastq.gz'
			rpgz=makeTrimmedFastQDir+'/'+sraIDs+'_2P.fastq.gz'
			rugz=makeTrimmedFastQDir+'/'+sraIDs+'_2U.fastq.gz'
			Trimmomatic_Command= "trimmomatic PE %s %s %s %s %s %s SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:%s:2:40:15 -threads %s" %(fastQ1, fastQ2, fpgz, fugz, rpgz, rugz, str(args.adapterPE), core)
			if args.verbose:
				print(Trimmomatic_Command)
			Trimmomatic_Command_run=subprocess.Popen(Trimmomatic_Command, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
			Trimmomatic_Command_run.wait()
	else: #SINGLE
		if sraIDs in prevTrim:
			pass
		else:
			fastQS=(glob.glob(makedataDir+'/'+sraIDs+fastQ_S))[0]
			fsgz=makeTrimmedFastQDir+'/'+sraIDs+'_1S.fastq.gz'
			Trimmomatic_Command= "trimmomatic SE %s %s SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:%s:2:40:15 -threads %s" %(fastQS, fsgz, str(args.adapterSE), core)
			if args.verbose:
				print(Trimmomatic_Command)
			Trimmomatic_Command_run=subprocess.Popen(Trimmomatic_Command, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
			Trimmomatic_Command_run.wait()


#FastQC for Trimmed data

prevTrimFastQC=set()
makeTrimFastQCDir= './Trimmed_Data_FastQC_report'
if not os.path.isdir(makeTrimFastQCDir):
	os.mkdir(makeTrimFastQCDir)
else:
	for items in (glob.glob(makeTrimFastQCDir+'/'+'*html')):
		if items:
			fileName=items.split('/')[-1].split('_')[0]
			prevTrimFastQC.add(fileName)
		else:
			fileName='Null'
			prevTrimFastQC.add(fileName)

TrimmedfastQC_com=[]


for sraIDs in sraSet:
	if sraSet[sraIDs][0]=='P': #PAIRED
		if sraIDs in prevTrimFastQC:
			pass
		else:
			trimfastQC_Command="fastqc %s -o %s" %(makeTrimmedFastQDir+'/'+sraIDs+'*P.*', makeTrimFastQCDir)
			TrimmedfastQC_com.append(trimfastQC_Command)
	else:
		if sraIDs in prevTrimFastQC:
			pass
		else:
			trimfastQC_Command="fastqc %s -o %s" %(makeTrimmedFastQDir+'/'+sraIDs+'*S.*', makeTrimFastQCDir)
			TrimmedfastQC_com.append(trimfastQC_Command)

print('\n>> FastQC Running for Trimmed Data .. ')

if TrimmedfastQC_com:
	thread_count = core # maximum parallel threads
	stopped = threading.Event()

	q = Queue()
	if args.verbose:
		print('-- Processing QC with FastQC for Trimmed Data: '+ (str(len(TrimmedfastQC_com))+ ' tasks in thread queue with '+ str(thread_count))+ ' thread limit')
	for item in TrimmedfastQC_com:
		q.put(item)
	for x in range(thread_count):
		t = threading.Thread(target=worker_func)
		# t.daemon = True #Enable to run threads as daemons
		t.start()
	q.join()       # block until all tasks are done
	stopped.set()

# multiQC run for Trimmed Data

makeTrimMultiQCInitialDir= './Trimmed_Data_multiQC_report'
if not os.path.isdir(makeTrimMultiQCInitialDir):
	os.mkdir(makeTrimMultiQCInitialDir)

trim_multiQC_command_init="multiqc -f %s -o %s" %(makeTrimFastQCDir, makeTrimMultiQCInitialDir)
print('\n>> MultiQC Running for Trimmed Data .. ')
trim_multiQC_command_run=subprocess.Popen(trim_multiQC_command_init, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
trim_multiQC_command_run.wait()

#Assembly Script

#for sraIDs in sraSet:
#	mapTrinity_Dir=sraIDs
#	if sraSet[sraIDs][0]=='P': #PAIRED
#		pgz=makeTrimmedFastQDir+'/'+sraIDs+'_*P.fastq.gz'
#		gunzip_Command="gunzip %s" %(pgz)
#		print(gunzip_Command)
#		gunzip_Command_run=subprocess.Popen(gunzip_Command, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
#		gunzip_Command_run.wait()
#	else:
#		pgz=makeTrimmedFastQDir+'/'+sraIDs+'_*S.fastq.gz'
#		gunzip_Command="gunzip %s" %(pgz)
#		print(gunzip_Command)
#		gunzip_Command_run=subprocess.Popen(gunzip_Command, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
#		gunzip_Command_run.wait()

print('\n'+'>> Assembly will start now... '+'\n\n')


mapSpades_Dir= './MappingSpades'
if not os.path.isdir(mapSpades_Dir):
	os.mkdir(mapSpades_Dir)

mapMegaHit_Dir='./MappingMegahit'
if not os.path.isdir(mapMegaHit_Dir):
	os.mkdir(mapMegaHit_Dir)

for sraIDs in sraSet:
	if sraSet[sraIDs][0]=='P': #PAIRED
		fpgz=makeTrimmedFastQDir+'/'+sraIDs+'_1P.fastq.gz'
		rpgz=makeTrimmedFastQDir+'/'+sraIDs+'_2P.fastq.gz'

		#trinity_Command="Trinity --seqType fq --left %s --right %s --CPU %s --output %s" %(fpgz, rpgz, core, mapTrinity_Dir)
		#print(trinity_Command)
		#trinity_Command_run=subprocess.Popen(trinity_Command, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
		#trinity_Command_run.wait()

		spades_Command="spades.py -1 %s -2 %s -o %s --rna -t %s"%(fpgz, rpgz, mapSpades_Dir+'/'+sraIDs, core)
		print(spades_Command)
		spades_Command_run=subprocess.Popen(spades_Command, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
		spades_Command_run.wait()

		megahit_Command="megahit -1 %s -2 %s -o %s --out-prefix %s -t %s"%(fpgz, rpgz, mapMegaHit_Dir+'/'+sraIDs, sraIDs, core)
		print(megahit_Command)
		megahit_Command_run=subprocess.Popen(megahit_Command, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
		megahit_Command_run.wait()

	else:
		fsgz=makeTrimmedFastQDir+'/'+sraIDs+'_1S.fastq.gz'

		#trinity_Command="Trinity --seqType fq --single %s --CPU %s --output %s" %(fsgz, core, mapTrinity_Dir)
		#print(trinity_Command)
		#trinity_Command_run=subprocess.Popen(trinity_Command, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
		#trinity_Command_run.wait()

		spades_Command="spades.py -s %s -o %s --rna -t %s"%(fsgz, mapSpades_Dir+'/'+sraIDs, core)
		print(spades_Command)
		spades_Command_run=subprocess.Popen(spades_Command, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
		spades_Command_run.wait()

		megahit_Command="megahit -r %s -o %s --out-prefix %s -t %s"%(fsgz, mapMegaHit_Dir+'/'+sraIDs, sraIDs, core)
		print(megahit_Command)
		megahit_Command_run=subprocess.Popen(megahit_Command, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
		megahit_Command_run.wait()



sys.exit()
