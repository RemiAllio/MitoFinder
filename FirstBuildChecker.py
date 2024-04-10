#!/usr/bin/python2.7
#Version: 1.0
#Author: Alex Schomaker - alexschomaker@ufrj.br
#LAMPADA - IBQM - UFRJ

'''
Copyright (c) 2014 Alex Schomaker Bastos - LAMPADA/UFRJ

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

from Bio import SeqIO, SearchIO
from Bio.Alphabet import generic_dna, generic_protein
from subprocess import Popen
import shlex, os
import tRNAscanChecker, circularizationCheck, geneChecker
from bisect import bisect_left

def checkResults(processName, pathToWork, sizeToLook, refSeqFile = None, cutoffValue = (2500,18000), blasteVal = 10.0,
                    blastHitSizePercentage = 0.625, usingSOAP = True, numberOfReadGroups = 1, buildCloroplast = False,
                    skipTrnaScan = False, circularSize = 40, circularOffSet = 220, cutoffEquality = 0.625, organismType = 2,
                    blastFolder = 'installed', noExtension = False, ignoreFirstBuildChecks = False, coveCutOff = 8,
                    buildBacteria = False, buildArchea = False, validContigs = 1):
	'''
	Do checks for main checker function...
	Checks for tRNAs, other features like rRNAs, genes, etc and returns a Assembly object, from tRNAscanChecker
	that holds this info.
	'''

	genBankReference = False
	
	#Let's set the default gene checking genbank file according to the organismType flag
	#and if a genbank reference was given in the -r flag, it will be changed
	module_dir = os.path.dirname(__file__)
	module_dir = os.path.abspath(module_dir)
	'''
		1. The Standard Code
		2. The Vertebrate Mitochondrial Code
		3. The Yeast Mitochondrial Code
		4. The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code
		5. The Invertebrate Mitochondrial Code
		6. The Ciliate, Dasycladacean and Hexamita Nuclear Code
		9. The Echinoderm and Flatworm Mitochondrial Code
		10. The Euplotid Nuclear Code
		11. The Bacterial, Archaeal and Plant Plastid Code
		12. The Alternative Yeast Nuclear Code
		13. The Ascidian Mitochondrial Code
		14. The Alternative Flatworm Mitochondrial Code
		16. Chlorophycean Mitochondrial Code
		21. Trematode Mitochondrial Code
		22. Scenedesmus obliquus Mitochondrial Code
		23. Thraustochytrium Mitochondrial Code
		24. Pterobranchia Mitochondrial Code
		25. Candidate Division SR1 and Gracilibacteria Code
	'''
	if buildCloroplast == True:
		refSeqFileForGenes = os.path.join(module_dir, 'references/cloroplast.gb')
	elif organismType == 11: #plant plastid or bacterial or archea
		refSeqFileForGenes = os.path.join(module_dir, 'references/magnolia.gb')
		if buildBacteria == True:
			refSeqFileForGenes = os.path.join(module_dir, 'references/bacteria.gb') #some bacterial DNA
		elif buildArchea == True:
			refSeqFileForGenes = os.path.join(module_dir, 'references/archea.gb') #some archea DNA
	elif organismType == 3: #yeast
		refSeqFileForGenes = os.path.join(module_dir, 'references/yeast.gb')
	elif organismType == 5: #invertebrate, insect
		refSeqFileForGenes = os.path.join(module_dir, 'references/beetle.gb')
	elif organismType == 6: #ciliate
		refSeqFileForGenes = os.path.join(module_dir, 'references/paramecium.gb')
	elif organismType != 1: #human, 2(vertebrate), also the default option
		refSeqFileForGenes = os.path.join(module_dir, 'references/human.gb')
	
	#we don't need a sequence reference for this check
	if refSeqFile != None: #if user gave a reference file not in fasta, let's consider it to look for features
		if refSeqFile[-6:] != '.fasta':
			genBankReference = True
			refSeq = SeqIO.read(refSeqFile, "genbank", generic_dna)
			refSeqFileForGenes = refSeqFile

	circularCheck = circularizationCheck.circularizationCheck("best_query.fasta", circularSize, circularOffSet, blastFolder) #returns a tuple with True or False and coordinates
	if circularCheck[0] == True:
		print('Circularization was found.')
	else:
		print('Circularization was not found.')

	if skipTrnaScan == False: #user did not want to skip tRNAscanChecker
		print("Checking for tRNAs...")
	else:
		print("Checking for genomic features...")
	checktRNA = tRNAscanChecker.tRNAscanCheck("best_query.fasta", circularCheck, skipTrnaScan, organismType, coveCutOff, buildBacteria,
						buildArchea)
	#add to the assembly object the number of contigs concatenated to create this super contig
	checktRNA.validContigs = validContigs
	'''
	let's check for it's features to see if everything was built,
	if a .gb reference was given we check against that, if not, we check against
	our own database inside references/ folder according to -o flag
	'''
	if ignoreFirstBuildChecks == False:
		checktRNA.checkFeatures = geneChecker.geneCheck(refSeqFileForGenes, "best_query.fasta", cutoffEquality, genBankReference, blastFolder,
                                                                organismType = organismType)
		presentFeatures = checktRNA.checkFeatures[0]
		importantFeatures = checktRNA.checkFeatures[1]
		splitFeatures = checktRNA.checkFeatures[2]
		print 'Features found: %s / %s' % (len(presentFeatures) - len(splitFeatures), len(importantFeatures))
	else:
		checktRNA.checkFeatures = ([],[])
		presentFeatures = checktRNA.checkFeatures[0]
		importantFeatures = checktRNA.checkFeatures[1]
		print 'Ignoring genomic checks, since --relaxed is on.'

	print "Assessing checks..."
	if (len(checktRNA) == 21 and buildCloroplast == False) or (len(checktRNA) == 37 and buildCloroplast == True):
		print 'All ' + str(len(checktRNA)) + ' tRNAs were found!'
		if circularCheck[0] == True:
			print 'And this build is circular.'
			if len(presentFeatures) == len(importantFeatures):
				print 'And all genes, tRNAs and rRNAs were found.'
				if skipTrnaScan == True:
					return True
			elif len(presentFeatures) >= len(importantFeatures):
				print 'And all genes, tRNAs and rRNAs were found, but some were split or duplicated.'
			else:
				print 'But not all genes, tRNAs and rRNAs were found, storing info and rebuilding...'
				print ''
				return checktRNA
			print ''
			#if all tRNAs were found, all features were found and it's circular, just return true and stop the recursive part
			return True
		else:
			print 'But circularization could not be found yet, storing this info...'
			if len(presentFeatures) == len(importantFeatures):
				print 'And all genes, tRNAs and rRNAs were found.'
			elif len(presentFeatures) >= len(importantFeatures):
				print 'And all genes, tRNAs and rRNAs were found, but some were split or duplicated.'
			else:
				print 'But not all genes, tRNAs and rRNAs were found, storing info...'
			print 'Rebuilding...'
			print ''
			return checktRNA
	else:
		print str(len(checktRNA)) + ' tRNAs were built.'
		if len(presentFeatures) == len(importantFeatures):
			print 'And all genes and rRNAs were found.'
		elif len(presentFeatures) >= len(importantFeatures):
			print 'And all genes, tRNAs and rRNAs were found, but some were split or duplicated.'
		else:
			print 'But not all genes and rRNAs were found, storing info...'
		print ''
		return checktRNA

def checkSoapOutput(processName, pathToWork, sizeToLook, refSeqFile = None, cutoffValue = (2500,18000), blasteVal = 10.0,
                    blastHitSizePercentage = 0.60, usingSOAP = True, numberOfReadGroups = 1, buildCloroplast = False,
                    skipTrnaScan = False, circularSize = 50, circularOffSet = 220, cutoffEquality = 0.60, organismType = 2,
                    blastFolder = 'installed', noExtension = False, ignoreFirstBuildChecks = False, coveCutOff = 8,
                    buildBacteria = False, buildArchea = False):
	'''
	DeNovo result checker. It tries to find out if a possible sequence was built by checking size and Blasting.
	'''
	#read SOAPdenovo-Trans scafSeq file and see if there is a possible target DNA built
	genBankReference = False
	if refSeqFile != None: #if user gave a reference file, let's consider its size as sizeToLook
		if refSeqFile[-6:] != '.fasta':
			genBankReference = True
			refSeq = SeqIO.read(refSeqFile, "genbank", generic_dna)
			#create fasta file for blast+ to work with
			blastFastaFile = refSeqFile + '.fasta'
			SeqIO.write(refSeq, blastFastaFile, 'fasta')
			refSeqFileForGenes = refSeqFile
			sizeToLook = len(refSeq.seq)
		else:
			refSeq = SeqIO.read(refSeqFile, "fasta", generic_dna)
			blastFastaFile = refSeqFile
			sizeToLook = len(refSeq)
	
	print ''
	print "Checking DeNovo results for possible hits!\nLogs being saved to " + pathToWork
	if usingSOAP == True:
		scafFile = pathToWork + processName + '.scafSeq'
	elif usingSOAP == 'Spades':
		scafFile = pathToWork + processName + '/scaffolds.fasta'
	else:
		if numberOfReadGroups > 1:
			scafFile = pathToWork + processName + '-denovo_assembly/' + processName + '-denovo_d_results/' + processName + '-denovo_out_AllStrains.unpadded.fasta'
		else:
			scafFile = pathToWork + processName + '-denovo_assembly/' + processName + '-denovo_d_results/' + processName + '-denovo_out.unpadded.fasta'

	possible_sequences = [] # Setup an empty list
	listOfValidResults = []
	minSizeToLook = cutoffValue[0]
	maxSizeToLook = cutoffValue[1]

	for record in SeqIO.parse(open(scafFile, "rU"), "fasta", generic_dna):
		if len(record.seq) >= minSizeToLook and len(record.seq) <= maxSizeToLook:
			# Add this record to our list
			possible_sequences.append(record)
 
	print "Found %i possible sequences due to size." % len(possible_sequences)
	
	possible_sequences = sorted(possible_sequences)
	output_handle = open("possible_hits.fasta", "w")
	SeqIO.write(possible_sequences, output_handle, "fasta")
	output_handle.close()
	
	#Blast part of checker!
	if len(possible_sequences) == 0:
		print "No possible sequences were found according to size check. Skipping the other checks."
		print ''
		return False
	elif refSeqFile != None:
		print "Formatting database for blast..."
		if blastFolder == 'installed':
			command = "formatdb -i " + blastFastaFile + " -p F" #need to formatdb refseq first
		else:
			command = blastFolder + "/bin/makeblastdb -in " + blastFastaFile + " -dbtype nucl" #need to formatdb refseq first

		args = shlex.split(command)
		formatDB = Popen(args, stdout=open(os.devnull, 'wb'))
		formatDB.wait()
		
		print "Running blast against refSeq to determine if a hit was built..."
		with open("possible_hits.blast.xml",'w') as blastResultFile:
			if blastFolder == 'installed':
				command = "blastall -p blastn -d " + refSeqFile + " -i possible_hits.fasta -e " + str(blasteVal) + " -m 7" #call BLAST with XML output
			else:
				command = blastFolder + "/bin/blastn -task blastn -db " + blastFastaFile + " -query possible_hits.fasta -evalue " + str(blasteVal) + " -outfmt 5 -num_threads 2" #call BLAST with XML output
			args = shlex.split(command)
			blastAll = Popen(args, stdout=blastResultFile)
			blastAll.wait()
		
		#SearchIO object handler and checker for best hit separation
		blastparse = SearchIO.parse('possible_hits.blast.xml', 'blast-xml') #get all queries
		dictOfBlastResults = {}
		for qresult in blastparse: #in each query, let's check for span size
			totalSpan = 0
			for hsp in qresult.hsps: #let's sum up all HSPs span sizes
				totalSpan += hsp.aln_span
			#if total blasted region is bigger or equal to minimum size to look times blasthitpercentage
			#add it to the dictionary of results
			if totalSpan >= minSizeToLook * blastHitSizePercentage:
				dictOfBlastResults[totalSpan] = qresult

		#sort the dicionary keys, which are the span sizes, in decrescent order
		orderedSpanSizes = sorted(dictOfBlastResults.keys(), reverse=True)

		#get the best result, the one which aligned more and use it as backbone for the rest of the process
		if len(orderedSpanSizes) > 0:
			qresult = dictOfBlastResults[orderedSpanSizes[0]]
			listOfValidResults.append(qresult)
			print("Best contig found: %s" % qresult.id)
		else:
			print 'Could not find a sequence that is good enough...\n'
			return False

		'''
		Down here we are going to try and increase the final sequence, by looking through the contigs blasted
		and seeing if they complement eachother. Doing this only if the best build found is lower than 92.5% of
		target size. Otherwise, just grab best hit and procceed.
		'''
		for record in SeqIO.parse(open("possible_hits.fasta", "rU"), "fasta", generic_dna): #check all possible hits
			if record.id == qresult.id: #found the best match
				output_handle = open("best_query.fasta", "w")
				final_Record = record
				
				if len(final_Record.seq) < sizeToLook * 0.925 and noExtension == False: #if lower than 92.5%, try to increase this sequence
					print 'Trying to find other contigs that match the target reference...\nSize before extension: ', len(final_Record.seq)
					totalSize = 0 #totalSize will be increased until we get as far as 1.05 * the reference sequence in size
					extraSeqsFound = '' #the final sequence is going to be built here
					#triesToExtend = 0 #limiting to 4 tries, so we don't overdo it

					#hash of start positions vs sequence and a second hash to see if it was in the -1 strand
					dictOfStarts = {}
					dictOfComplements = {}
					dictOfEndings = {}
					dictOfQueryStarts = {}
					dictOfQueryEndings = {}

					#actual search down here
					blastparse2 = SearchIO.parse('possible_hits.blast.xml', 'blast-xml') #open xml parse
					for extendingResult in blastparse2: #in each query, let's check for span size
						addToDict = False #create this as false, if there is an alignment of at least 200 in size make this true
						if extendingResult.seq_len >= 150:
							for hsp_2 in extendingResult.hsps: #let's sum up all HSPs span sizes
								if hsp_2.aln_span >= 150:
									#to add to the dictOfComplements, in case we needed to reverse this seq
									reverseComplement = hsp_2.hit_range[0] > hsp_2.hit_range[1] or hsp_2.hit_frame == -1
									lowestStart = min(hsp_2.hit_range[0],hsp_2.hit_range[1])
									highestEnding = max(hsp_2.hit_range[0],hsp_2.hit_range[1])
									queryStart = min(hsp_2.query_range[0],hsp_2.query_range[1])
									queryEnding = max(hsp_2.query_range[0],hsp_2.query_range[1])
									addToDict = True
									totalSize += extendingResult.seq_len
									break

						if addToDict == True:
							#add this result to the hash, since it met the criteria
							if lowestStart not in dictOfStarts:
								dictOfStarts[lowestStart] = extendingResult
								dictOfEndings[lowestStart] = highestEnding
								dictOfComplements[lowestStart] = reverseComplement
								dictOfQueryStarts[lowestStart] = queryStart
								dictOfQueryEndings[lowestStart] = queryEnding
							else:
								if extendingResult.seq_len > dictOfStarts[lowestStart].seq_len:
									dictOfStarts[lowestStart] = extendingResult
									dictOfEndings[lowestStart] = highestEnding
									dictOfComplements[lowestStart] = reverseComplement
									dictOfQueryStarts[lowestStart] = queryStart
									dictOfQueryEndings[lowestStart] = queryEnding
						if totalSize > sizeToLook * 1.05:
							break

					'''
					Down here we clean up the blast results to make sure we don't erroneously extend the contigs.
					It was needed because some sequences blasted inside another one, but, due to assembly errors
					were bigger and had duplicated regions and were inserted improperly.
					'''
					toDeleteList = []
					for myStart in dictOfStarts:
						myEnding = dictOfEndings[myStart]
						for otherStart in dictOfStarts:
							if otherStart != myStart:
								otherEnding = dictOfEndings[otherStart]
								if otherStart <= myStart and otherEnding >= myEnding:
									toDeleteList.append(myStart)
									break
					
					for deleteMe in toDeleteList:
						del dictOfStarts[deleteMe]
						del dictOfEndings[deleteMe]
						del dictOfComplements[deleteMe]
						del dictOfQueryStarts[deleteMe]
						del dictOfQueryEndings[deleteMe]

					'''
					Sorting the list of results to make sure we increase from start to end.
					'''
					listOfValidResults = sorted(dictOfStarts.keys()) #list of found contigs sorted based on start positions
					numberOfReverseComplementsAdded = 0 #if this is > 0 then reverse final result, otherwise, don't

					'''
					Down here is where the increasing of the final De Novo sequence will happen.
					Looks for start position (based on blast against reference), and inserts it into the final
					sequence based on that.
					If it is insided an already covered region, ignore it.
					'''
					lastInsertEnding = 0
					for n in xrange(len(listOfValidResults)):
						print 'Found contig/scaffold with id: ', dictOfStarts[listOfValidResults[n]].id
						for record2 in SeqIO.parse(open("possible_hits.fasta", "rU"), "fasta", generic_dna):
							startVal = listOfValidResults[n]
							blastResult = dictOfStarts[startVal]
							if record2.id == blastResult.id:
								reverseComplement = dictOfComplements[startVal]
								queryStart = dictOfQueryStarts[startVal]
								queryEnd = dictOfQueryEndings[startVal]
								if reverseComplement == True:
									numberOfReverseComplementsAdded += 1
									seqToAppend = record2.seq.reverse_complement()[queryStart:queryEnd]
								else:
									numberOfReverseComplementsAdded -= 1
									seqToAppend = record2.seq[queryStart:queryEnd]
								if lastInsertEnding == 0:
									extraSeqsFound += seqToAppend
								elif startVal > lastInsertEnding:
									extraSeqsFound += 'n' * max(20,(startVal - lastInsertEnding)) + seqToAppend
								elif startVal <= lastInsertEnding:
									seqToAppend = seqToAppend[(-1) * (dictOfEndings[startVal] - lastInsertEnding):]
									newSequence = extraSeqsFound + seqToAppend
									if len(newSequence) > len(extraSeqsFound):
										extraSeqsFound = newSequence

								lastInsertEnding = dictOfEndings[startVal]

					#finished increasing, time to end it and report
					if len(extraSeqsFound) - extraSeqsFound.lower().count('n') > len(final_Record.seq):
						final_Record.seq = 'n'*20 + extraSeqsFound + 'n'*20
						#if numberOfReverseComplementsAdded > 0: #since we had to reverse many sequences to properly inser it
											#in the final sequence, reverse again to turn it back to original
							#final_Record.seq = final_Record.seq.reverse_complement()
						print 'Size after extension: ', len(final_Record.seq)
						print 'Size after extension (without Ns): ', len(final_Record.seq) - final_Record.seq.lower().count('n')
					else:
						print 'Extension was not possible, keeping original sequence.\n'
				
				SeqIO.write(final_Record, output_handle, "fasta")
				output_handle.close()
				'''
				#create a fasta file with everything but the best hit, needed later if edges are missing
				sequencesExceptBest = []
				for sequenceFound in SeqIO.parse(open(scafFile, "rU"), "fasta", generic_dna):
					if sequenceFound.id != qresult.id:
						sequencesExceptBest.append(sequenceFound)
				outputSeqs = open("all_hits_except_best.fasta", "w")
				SeqIO.write(sequencesExceptBest, outputSeqs, "fasta")
				outputSeqs.close()
				'''

				#Checking if there is circularization...
				print "A possible hit was found. Going to next step..."
				print 'Checking for circularization...'
				#call this function to check for circularization, genomic features and tRNAs
				return checkResults(processName, pathToWork, sizeToLook, refSeqFile, cutoffValue, blasteVal, blastHitSizePercentage, 
                                    usingSOAP, numberOfReadGroups, buildCloroplast, skipTrnaScan, circularSize, circularOffSet, 
                                    cutoffEquality, organismType, blastFolder, noExtension, ignoreFirstBuildChecks, coveCutOff,
						            buildBacteria, buildArchea, len(listOfValidResults))
	else: #user didn't provide a reference sequence to blast against, let's just consider the size then and move on!
		output_handle = open("best_query.fasta", "w")
		'''
		Do a bisection search to look for the contig with the closest size to target size, since no
		reference was given.
		'''
		pos = bisect_left(map(len,possible_sequences), sizeToLook)
		if pos == 0:
			SeqIO.write(possible_sequences[0], output_handle, "fasta")
		elif pos == len(possible_sequences):
			SeqIO.write(possible_sequences[-1], output_handle, "fasta")
		else:
			before = possible_sequences[pos - 1]
			after = possible_sequences[pos]
			if len(after) - sizeToLook < sizeToLook - len(before):
				SeqIO.write(after, output_handle, "fasta")
			else:
				SeqIO.write(before, output_handle, "fasta")
		output_handle.close()

		print "Possible hits were found. Grabbing sequence closest to optimum size... Going to next step..."
		print ''
		#check for circularization, tRNAs and genomic features
		#return True, False or Assembly class holding this build's information
		return checkResults(processName, pathToWork, sizeToLook, refSeqFile, cutoffValue, blasteVal, blastHitSizePercentage, usingSOAP,
                            numberOfReadGroups, buildCloroplast, skipTrnaScan, circularSize, circularOffSet, cutoffEquality, organismType,
                            blastFolder, noExtension, ignoreFirstBuildChecks, coveCutOff, buildBacteria, buildArchea, len(listOfValidResults))
