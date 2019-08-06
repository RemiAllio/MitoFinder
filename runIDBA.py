#!/usr/bin/python
#Version: 1.0.2
#Authors: Allio Remi & Schomaker-Bastos Alex
#ISEM - CNRS - LAMPADA - IBQM - UFRJ

'''
Copyright (c) 2019 Remi Allio - ISEM/CNRS & Alex Schomaker-Bastos - LAMPADA/UFRJ

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
'''

from subprocess import Popen
import gzip
import subprocess
import time
import shlex, os, shutil, FirstBuildChecker
from shutil import copyfile


def runIDBA(processName = 'teste', shortestContig = 100, inputFile = 'teste.input', processorsToUse = 4,
					idbaFolder= 'installed', refSeqFile = None, organismType = 2,
					blastFolder = 'installed', logfile='logfile'):

	pathToIdba = idbaFolder
	bestBuild = None

	print 'Starting Assembly step with IDBA-UD '
	logfile=open(logfile,"a")
	logfile.write('Starting Assembly phase with IDBA-UD '+"\n")
	
	pathToWork = os.getcwd()+"/"
	print 'Result files will be saved here: '
	print pathToWork+processName+"_idba/"
	logfile.write('Result files will be saved here: '+"\n"+pathToWork+processName+"_idba/"+"\n")
		
	#copy input file to secondary folder to keep it organized
	if "/" in inputFile:
		destFile = pathToWork + str(inputFile.split("/")[-1])
	else:
		destFile = pathToWork + "/" + inputFile
	shutil.copyfile(inputFile, destFile)

	#####################################
	######	    Run IDBA-UD!!!     ######
	#####################################
		
	out=processName+"_idba"
	idba="yes"
	if os.path.isdir(out):
		print "\n####################################"
		print "\n WARNING : "+pathToWork+out+" already exists." 
		print "Mitofinder will skip idba step"
		print "\nIf you want to run idba again, kill the mitofinder process, remove or rename the idba result folder, and restart mitofinder\n"
		print "#####################################\n"
		logfile.write("\n####################################"+"\n"+"\n WARNING : "+pathToWork+out+" already exists." +"\n"+"Mitofinder will skip idba step"+"\n"+"\nIf you want to run idba again, kill the mitofinder process, remove or rename the idba result folder, and restart mitofinder\n"+"\n"+"#####################################\n"+"\n")
		time.sleep(2)
		idba="no"
		exit()
	#create IDBA logfile:
	
	with open(inputFile,'r') as InputFile:
		with open(pathToWork + 'idba.log','w') as idbaLogFile:
			for line in InputFile:
					if '#' != line[0] and line != '\n':
						configPart = line.lower().replace('\n','').replace(' ','').split('=')[0]
						if configPart == 'type':
							t=line.replace('\n','').replace(' ','').split('=')[-1]
						if configPart == 'q1':
							read1=line.replace('\n','').replace(' ','').split('=')[-1]
							if not "/" in read1:
								read1="../"+read1
							if read1[-3:] == ".gz":
								with gzip.open(read1, 'rb') as f_in:
									with open(read1[0:-3], 'wb') as f_out:
										shutil.copyfileobj(f_in, f_out)
								read1=read1[0:-3]
						if configPart == 'q2':
							read2=line.replace('\n','').replace(' ','').split('=')[-1]
							if not "/" in read2:
								read2="../"+read2
							if read2[-3:] == ".gz":
								with gzip.open(read2, 'rb') as f_in:
									with open(read2[0:-3], 'wb') as f_out:
										shutil.copyfileobj(f_in, f_out)
								read2=read2[0:-3]

	if idba == "yes":
		with open(pathToWork + 'idba.log','a') as idbaLogFile:
			if t == "PE":
				print "Paired-end"
				logfile.write("Paired-end"+"\n")
				read=processName+"_idba_read.fasta"
				command = '%sfq2fa --merge --filter %s %s %s' %(pathToIdba, read1, read2, read)
				print "Preparing data for IDBA-UD assembly" 
				logfile.write("Preparing data for IDBA-UD assembly"+"\n")
				fq2fa = Popen(command, stdout=idbaLogFile, stderr=idbaLogFile, shell=True)
				fq2fa.wait()
				if not os.path.isfile(pathToWork+"/"+read) == True:
					print "\n ERROR: IDBA-UD didn't run well"
					print "Please check log file : "+ pathToWork + 'idba.log'
					logfile.write("\n ERROR: IDBA-UD didn't run well"+"\n"+"Please check log file : "+ pathToWork + 'idba.log'+"\n")
					exit()
				command = '%sidba -r %s -o %s --num_threads %s' %(pathToIdba, read, out, processorsToUse)
				print "Running assembly" 
				idba = Popen(command, stdout=idbaLogFile, stderr=idbaLogFile, shell=True)
				idba.wait()
				#copyfile(pathToWork+"/"+out+"/contig.fa", pathToWork+"/"+processName+".scafSeq")
				
			if t == "SE":
				print "Single-end"
				logfile.write("Single-end"+"\n")
				read=processName+"_idba_read.fasta"
				command = '%sfq2fa --filter %s %s' %(pathToIdba, read1, read)
				print "Preparing data for IDBA-UD assembly" 
				logfile.write("Preparing data for IDBA-UD assembly"+"\n")
				fq2fa = Popen(command, stdout=idbaLogFile, stderr=idbaLogFile, shell=True)
				fq2fa.wait()
				if not os.path.isfile(pathToWork+"/"+read) == True:
					print "\n ERROR : IDBA-UD didn't run well"
					print "Please check log file : "+ pathToWork + 'idba.log'
					logfile.write("\n ERROR: IDBA-UD didn't run well"+"\n"+"Please check log file : "+ pathToWork + 'idba.log'+"\n")
					exit()
				command = '%sidba -r %s -o %s --num_threads %s' %(pathToIdba, read, out, processorsToUse)
				print "Running assembly"
				logfile.write("Running assembly"+"\n") 
				idba = Popen(command, stdout=idbaLogFile, stderr=idbaLogFile, shell=True)
				idba.wait()
				#copyfile(pathToWork+"/"+out+"/contig.fa", pathToWork+"/"+processName+".scafSeq")
	
	with gzip.open(read+'.gz', 'wb') as f:
		f.write(read)
	
	logfile.close()








































