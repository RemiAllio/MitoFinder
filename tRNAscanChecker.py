#!/usr/bin/python
#Version: 1.4
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

from Bio import SeqIO
from Bio.Alphabet import generic_dna, generic_protein, IUPAC
from Bio.Data import CodonTable
from subprocess import Popen
import shlex, sys, os

class Assembly():
	'''
	Class to hold the sequence from best_query and it's tRNAscan-SE results, as well as features.
	'''
	def __init__(self, resultFile, tRNAscanResultFile, hasCircularized, tRNAscan, organismType):
		if resultFile[-6:] != '.fasta':
			self.refSeq = SeqIO.read(resultFile, "genbank", generic_dna)
		else:
			self.refSeq = SeqIO.read(resultFile, "fasta", generic_dna)
		self.tRNAs = []
		self.checkFeatures = None
		self.validContigs = 1
		checkedtRNAs = []
		self.hasCircularized = hasCircularized
		startCheck = False
		if tRNAscanResultFile != None:
			with open(tRNAscanResultFile, 'r') as tRNAscanFile:
				if tRNAscanResultFile[-6:] == '.mitfi':
					module_dir = os.path.dirname(__file__)
					module_dir = os.path.abspath(module_dir)
					module_dir = os.path.join(module_dir, "mitfi/")
					dico_tRNA={}
					for line in open(module_dir+"Codon-Amino_Acid_Abbreviations.csv"):
						dico_tRNA[line.rstrip().split(";")[3]]=line.split(";")[2]
					c=1
					for line in tRNAscanFile:
						if line[:7] == '#header': 
							'''
							look for the start of tRNA results as in:
							#header	start	stop	evalue	AC	AA	model	strand
							'''
							startCheck = True
						elif startCheck == True and line[0] != "#": #we found the start, so let's start saving this info!
							tRNAinfos = line.split()
							#sequence name is tRNAinfos[0], but isn't really needed
							tRNAnumber = c
							c+=1
							tRNAstrand=tRNAinfos[8]
							if tRNAstrand == "+":
								tRNAcoordinates = (int(tRNAinfos[1]),int(tRNAinfos[2]))
							elif tRNAstrand == "-":
								tRNAcoordinates = (int(tRNAinfos[2]),int(tRNAinfos[1]))
							tRNAtype = str(tRNAinfos[6])
							tRNAtype = dico_tRNA.get(tRNAtype[0])
							if tRNAtype in checkedtRNAs:
								tRNAtype += str(checkedtRNAs.count(tRNAtype) + 1)
							checkedtRNAs.append(tRNAtype)
							tRNAcodon = tRNAinfos[5]
							tRNAintronBegin = 0
							tRNAintronEnd = 0
							if tRNAintronBegin > 0:
								print '	WARNING: %s found with an intron. Strange result, check .trnascan file' % tRNAtype
							tRNAscore = float(tRNAinfos[3])
							thisRNA = self.tRNA(self, tRNAnumber, tRNAcoordinates, tRNAtype, tRNAcodon, tRNAscore, tRNAintronBegin, tRNAintronEnd)
							self.tRNAs.append(thisRNA)
				elif tRNAscanResultFile[-9:] == '.trnascan':
					for line in tRNAscanFile:
						if line[0] == '-': 
							'''
							look for the start of tRNA results as in:
							tRNA #	Begin	End  	Type	Codon	Begin	End	Score
							'''
							startCheck = True
						elif startCheck == True: #we found the start, so let's start saving this info!
							tRNAinfos = line.split()
							#sequence name is tRNAinfos[0], but isn't really needed
							tRNAnumber = int(tRNAinfos[1])
							tRNAcoordinates = (int(tRNAinfos[2]),int(tRNAinfos[3]))
							tRNAtype = tRNAinfos[4]
							if tRNAtype in checkedtRNAs:
								tRNAtype += str(checkedtRNAs.count(tRNAtype) + 1)
							checkedtRNAs.append(tRNAinfos[4])
							tRNAcodon = tRNAinfos[5]
							tRNAintronBegin = int(tRNAinfos[6])
							tRNAintronEnd = int(tRNAinfos[7])
							if tRNAintronBegin > 0:
								print '	WARNING: %s found with an intron. Strange result, check .trnascan file' % tRNAtype
							tRNAscore = float(tRNAinfos[8])
							thisRNA = self.tRNA(self, tRNAnumber, tRNAcoordinates, tRNAtype, tRNAcodon, tRNAscore, tRNAintronBegin, tRNAintronEnd)
							self.tRNAs.append(thisRNA)
				elif tRNAscanResultFile[-6:] == '.arwen':
					for line in tRNAscanFile:
						if "genes found" in line: 
							'''
							look for the start of tRNA results as in:
							tRNA #	Begin	End  	Type	Codon	Begin	End	Score
							
							arwen
							# mtRNA-Type [Begin,End] N (codon)
												
							
							'''
							startCheck = True
						elif startCheck == True: #we found the start, so let's start saving this info!
							tRNAinfos = line.split()
							#sequence name is tRNAinfos[0], but isn't really needed
							tRNAnumber = int(line.split(" ")[0])
							if "c[" in line :
								tRNAcoordinates = (line.split("]")[0].split(",")[1],line.split("[")[1].split(",")[0])
							else:
								tRNAcoordinates = (line.split("[")[1].split(",")[0],line.split("]")[0].split(",")[1])
							tRNAtype = line.split("-")[1].split(" ")[0]
							if tRNAtype in checkedtRNAs:
								tRNAtype += str(checkedtRNAs.count(tRNAtype) + 1)
							checkedtRNAs.append(tRNAtype)
							tRNAcodon = line.split("(")[1].split(")")[0].upper()
							tRNAintronBegin = 0
							tRNAintronEnd = 0
							tRNAscore = 50.0
							thisRNA = self.tRNA(self, tRNAnumber, tRNAcoordinates, tRNAtype, tRNAcodon, tRNAscore, tRNAintronBegin, tRNAintronEnd)
							self.tRNAs.append(thisRNA)
		
	def __len__(self):
		return len(self.tRNAs)

	def isCircular(self):
		return self.hasCircularized
		
	class tRNA():
		'''
		Class to hold the tRNAscan-SE info for each tRNA in Assembly()
		'''
		def __init__(self, motherSeq, tRNAnumber, tRNAcoordinates, tRNAtype, tRNAcodon, tRNAscore, tRNAintronBegin, tRNAintronEnd):
			self.tRNAnumber = tRNAnumber
			reverseStart = False
			if int(tRNAcoordinates[0]) > int(tRNAcoordinates[1]):
				reverseStart = True
			tRNAstart = min(int(tRNAcoordinates[0]),int(tRNAcoordinates[1]))
			tRNAend = max(int(tRNAcoordinates[0]),int(tRNAcoordinates[1]))
			tRNAstart -= 1
			tRNAstart = max(0,tRNAstart)
			if reverseStart is False:
				tRNAcoordinates = [tRNAstart,tRNAend]
			else:
				tRNAcoordinates = [tRNAend,tRNAstart]
			self.tRNAcoordinates = tRNAcoordinates
			self.tRNAtype = tRNAtype
			self.tRNAcodon = tRNAcodon
			self.tRNAscore = tRNAscore
			self.tRNAintronBegin = tRNAintronBegin
			self.tRNAintronEnd = tRNAintronEnd
			self.motherSeq = motherSeq
			
		def __len__(self):
			return self.coordinates()[1] - self.coordinates()[0]
			
		def coordinates(self):
			return self.tRNAcoordinates
		
		def number(self):
			return self.tRNAnumber
		
		def typeOfRna(self):
			return self.tRNAtype
		
		def codon(self):
			return self.tRNAcodon
		
		def score(self):
			return self.tRNAscore

def tRNAscanCheck(resultFile = None, hasCircularized = False, skipTRNA = False, organismType = 2, coveCutOff = 7,
                  buildBacteria = False, buildArchea = False, tRNAscan="mitfi"):
	'''
	Use tRNAscan-SE to look for tRNAs and hold it's positions and scores in the tRNA Class.
	'''
	if skipTRNA == False:
		module_dir = os.path.dirname(__file__)
		module_dir = os.path.abspath(module_dir)
		scanInput = resultFile
		if tRNAscan == "mitfi":
			module_dir = os.path.join(module_dir, "mitfi/")
			outputName = resultFile[0:-6] + ".mitfi"
		elif tRNAscan == "trnascan":
			module_dir = os.path.join(module_dir, "trnascanSE/tRNAscan-SE-2.0/")
			outputName = resultFile[0:-6] + ".trnascan"
		elif tRNAscan == "arwen":
			module_dir = os.path.join(module_dir, "arwen/")
			outputName = resultFile[0:-6] + ".arwen"
		cfg_dir = os.path.dirname(__file__)
		cfg_full_path = os.path.join(cfg_dir, 'Mitofinder.config')

		with open(cfg_full_path,'r') as configFile:
			#grabbing the tRNAscan-SE folder from the config file...
			for line in configFile:
				if '#' != line[0] and line != '\n':
					configPart = line.lower().replace('\n','').replace(' ','').split('=')[0]
					if configPart == 'mitfifolder':
						MitFiFolder = line.replace('\n','').replace(' ','').split('=')[-1]
						if MitFiFolder.lower() == 'default':
							MitFiFolder = module_dir
					elif configPart == 'arwenfolder':
						ArwenFolder = line.replace('\n','').replace(' ','').split('=')[-1]
						if ArwenFolder.lower() == 'default':
							ArwenFolder = module_dir
					elif configPart  == 'trnascanfolder':
						tRNAscanFolder = line.replace('\n','').replace(' ','').split('=')[-1]
						if tRNAscanFolder.lower() == 'default':
							tRNAscanFolder = module_dir
		#adding the Arwen folder to these environments in the OS to avoid errors being thrown by Arwen
		if tRNAscan == "mitfi":	
			try:
				os.environ["PATH"] += os.pathsep + MitFiFolder
				os.environ["PERL5LIB"] += os.pathsep + MitFiFolder
			except KeyError:
				pass
			except:
				print 'MitFi path is not correct. Change it in Mitofinder.config'
		elif tRNAscan == "trnascan":	
			try:
				os.environ["PATH"] += os.pathsep + tRNAscanFolder
				os.environ["PERL5LIB"] += os.pathsep + tRNAscanFolder
			except KeyError:
				pass
			except:
				print 'tRNAscan path is not correct. Change it in Mitofinder.config'
		elif tRNAscan == "arwen":	
			try:
				os.environ["PATH"] += os.pathsep + ArwenFolder
				os.environ["PERL5LIB"] += os.pathsep + ArwenFolder
			except KeyError:
				pass
			except:
				print 'Arwen path is not correct. Change it in Mitofinder.config'
			
		if os.path.exists(outputName): #remove result file if it already exists so that tRNAscan doesn't throw another error
			os.remove(outputName)

		if tRNAscan == "mitfi":
			command = "ln -s "+MitFiFolder+"cmsearch ./" 
			args = shlex.split(command)
			tRNAscanRun = Popen(args, cwd=MitFiFolder)
			try:
				with open(outputName,"w") as tRNAscanLog:
					if MitFiFolder.lower() == 'installed':
						command = "java -jar mitfi.jar -code "+ str(organismType) + " " + scanInput 
						args = shlex.split(command)
						tRNAscanRun = Popen(args, stdout=tRNAscanLog, stderr=tRNAscanLog)
					else:
						command = "java -jar "+ MitFiFolder +"mitfi.jar -code "+ str(organismType) + " " + scanInput 
						args = shlex.split(command)
						tRNAscanRun = Popen(args, cwd=MitFiFolder, stdout=tRNAscanLog, stderr=tRNAscanLog)
					tRNAscanRun.wait()
	
				thisSequenceResult = Assembly(resultFile, outputName, hasCircularized, tRNAscan, organismType)
				return thisSequenceResult
			except:
				print ''
				print "MitFi failed."
				print ''
				return False

		if tRNAscan == "trnascan":
			#down here, we check organism type to make tRNAscan-SE run in appropriate mode
			if buildBacteria == True:
				organismFlag = '-B ' #bacterial
			elif buildArchea == True:
				organismFlag = '-A ' #archea
			elif organismType == 1:
				organismFlag = '' #leave it empty to use eukariotyc search
			else:
				organismFlag = '-O '
	
			#use different genetic code?
			if organismType == 2: #vertebrate
				geneticCode = '-g lib/gcode/gcode.vertmito '
			elif organismType == 3: #yeast
				geneticCode = '-g lib/gcode/gcode.ystmito '
			elif organismType == 4: #mold, protozoan
				geneticCode = '-g lib/gcode/gcode.othmito '
			elif organismType == 5: #invertebrate
				geneticCode = '-g lib/gcode/gcode.invmito '
			elif organismType == 6: #ciliate
				geneticCode = '-g lib/gcode/gcode.cilnuc '
			elif organismType == 9: #echinodermata
				geneticCode = '-g lib/gcode/gcode.echdmito '
			else:
				geneticCode = ''
			
			try:
				with open("tRNAscan-SE.log","w") as tRNAscanLog:
					if tRNAscanFolder.lower() == 'installed':
						command = "tRNAscan-SE -X " + str(coveCutOff) + ' ' + geneticCode + organismFlag + "-o " + outputName + " " + scanInput
						args = shlex.split(command)
						tRNAscanRun = Popen(args, stdout=tRNAscanLog, stderr=tRNAscanLog)
					else:
						command = "tRNAscan-SE -X " + str(coveCutOff) + ' ' + geneticCode + organismFlag + "-o " + outputName + " " + scanInput
						args = shlex.split(command)
						tRNAscanRun = Popen(args, cwd=tRNAscanFolder, stdout=tRNAscanLog, stderr=tRNAscanLog)
					tRNAscanRun.wait()
	
				thisSequenceResult = Assembly(resultFile, outputName, hasCircularized, tRNAscan, organismType)
				return thisSequenceResult
			except:
				print ''
				print "tRNAscan-SE failed. Check it's logs for more details."
				print ''
				return False
		
		elif tRNAscan == "arwen":
			#down here, we check organism type to make ARWEN run in appropriate mode
			if buildBacteria == True:
				geneticCode = '-gcbact ' #bacterial
			elif buildArchea == True:
				geneticCode = '-gcbact ' #archea
			elif organismType == 1:
				geneticCode = '-gcstd ' #standard
			elif organismType == 2: #vertebrate
				geneticCode = '-gcvert '
			elif organismType == 3: #yeast
				geneticCode = '-gcyeast '
			elif organismType == 4: #mold, protozoan
				geneticCode = '-gcprot '
			elif organismType == 5: #invertebrate
				geneticCode = '-gcinvert '
			elif organismType == 6: #ciliate
				geneticCode = '-gcciliate '
			elif organismType == 9: #echinodermata
				geneticCode = '-gcflatworm '
			else:
				geneticCode = ''
			
			try:
				with open("ARWEN.log","w") as tRNAscanLog:
					if ArwenFolder.lower() == 'installed':
						command = "arwen " + ' ' + geneticCode + "-o "\
							+ outputName + " -w " + scanInput
						args = shlex.split(command)
						tRNAscanRun = Popen(args, stdout=tRNAscanLog, stderr=tRNAscanLog)
					else:
						command = "arwen " + ' ' + geneticCode + "-o "\
							+ outputName + " -w " + scanInput
						print command
						args = shlex.split(command)
						tRNAscanRun = Popen(args, cwd=ArwenFolder, stdout=tRNAscanLog, stderr=tRNAscanLog)
					tRNAscanRun.wait()
	
				thisSequenceResult = Assembly(resultFile, outputName, hasCircularized, tRNAscan, organismType)
				return thisSequenceResult
			except:
				print ''
				print "Arwen failed. Check log files for more details."
				print ''
				return False
	else: #if tRNA = True
		return Assembly(resultFile, None, hasCircularized)
		
def tRNAconvert(inputTRNA):
	'''
	Gets as input: a list of tRNAs in tRNA-Phe format or TRNF
	Returns: a list or a string (if input was only a string of tRNA)
		with the format tRNA-Phe
	DOES NOT ALTER THE FINAL NAME, just the checks
	'''
	returnString = False
	if type(inputTRNA) != list:
		inputTRNA = [str(inputTRNA)]
		returnString = True
	outputList = []
	for convertTRNA in inputTRNA:
		dictOftRNAs = {'trnf':'trna-phe', 'trnv':'trna-val',
				'trnl':'trna-leu', 'trni':'trna-ile',
				'trnq':'trna-gln', 'trnm':'trna-met',
				'trnw':'trna-trp', 'trna':'trna-ala',
				'trnn':'trna-asn', 'trnc':'trna-cys',
				'trny':'trna-tyr', 'trns':'trna-ser',
				'trnd':'trna-asp', 'trnk':'trna-lys',
				'trng':'trna-gly', 'trnr':'trna-arg',
				'trnh':'trna-his', 'trne':'trna-glu',
				'trnt':'trna-thr', 'trnp':'trna-pro'}
				
		thisTRNA = convertTRNA.lower()
		thisPartialTRNA = thisTRNA[0:4]
		if len(thisTRNA) <= 6:
			if thisPartialTRNA in dictOftRNAs:
				if len(thisTRNA) == 4:
					finalName = dictOftRNAs[thisTRNA]
				else:
					if thisTRNA[-1] == '1':
						finalName = dictOftRNAs[thisPartialTRNA]
					else:
						finalName = dictOftRNAs[thisPartialTRNA] + thisTRNA[-1]
		else:
			finalName = thisTRNA

		outputList.append(finalName)

	if returnString == False:
		return outputList
	else:
		return outputList[0]

def prettyRNAName(tRNAName):
	'''
	Just a cosmetic function to properly format tRNA results that come in lowercase to the
	tRNA-Phe format.
	'''
	return tRNAName.replace('trna','tRNA')[0:5] + tRNAName[tRNAName.index('-') + 1:].capitalize()
