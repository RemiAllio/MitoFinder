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

from Bio import SeqIO, SeqFeature
from Bio.Alphabet import generic_dna, generic_protein, IUPAC
from Bio.Seq import Seq
from Bio.Data import CodonTable
from decimal import Decimal

def genbankOutput(resultGbFile, resultFile, listOfFeaturesToOutput, buildCloroplast = False, dLoopSize = 800):
	'''
	Creates a genbank file based on a fasta file given (resultfile) and a list of features that the genbank
	file should present (listoffeaturestooutput)
	'''
	#creating the genbank file, not annotated, to be opened afterwards and have the features inserted
	with open(resultGbFile, "w") as outputResult:
		finalResults = SeqIO.read(open(resultFile, 'rU'), "fasta", generic_dna)
		finalResults.seq = finalResults.seq.upper()
		finalResults.name = finalResults.name[0:10] + '_draft'
		finalResults.id = finalResults.name[0:10] + '_draft'
		if len(finalResults.name) > 16: #has to 16 characters long at max, or else genbank file throws error
			finalResults.name = finalResults.name[0:16]
			finalResults.id = finalResults.id[0:16]
		count = SeqIO.write(finalResults, outputResult, "genbank")

	with open(resultGbFile, "rU") as outputResult: #opening the output file, this time to insert the features
		finalResults = SeqIO.read(outputResult, "genbank", generic_dna)
		#lastFeatureAlignment = None
		dLoopFound = False
		for thisFeatureAlignment in listOfFeaturesToOutput:
			# 1. Define a feature type as a text string
			main_feature_qualifiers = {} #create qualifiers dict where the name will be stored

			if 'trn' in thisFeatureAlignment.seq2.lower() or 'rrn' in thisFeatureAlignment.seq2.lower() \
			or 'ribosomal' in thisFeatureAlignment.seq2.lower() or 'rnr' in thisFeatureAlignment.seq2.lower():
				main_feature_qualifiers['product'] = thisFeatureAlignment.seq2
				if 'trn' in thisFeatureAlignment.seq2.lower():
					main_feature_type = "tRNA"
				else:
					main_feature_type = "rRNA"
			else:
				main_feature_qualifiers['gene'] = thisFeatureAlignment.seq2
				main_feature_type = "gene"
				
			main_start_pos = SeqFeature.ExactPosition(thisFeatureAlignment.startBase)
			main_end_pos = SeqFeature.ExactPosition(thisFeatureAlignment.endBase)
			
			if main_feature_type == "gene":
				codonDiff = ((main_end_pos - (main_start_pos + 1)) % 3)
				if codonDiff == 2:
					main_end_pos += 1
				elif codonDiff == 1:
					main_end_pos -= 1

			# 2. Use the locations do define a FeatureLocation
			if thisFeatureAlignment.frame < 0:
				strandToOutput = -1
			else:
				strandToOutput = 1
			main_feature_location = SeqFeature.FeatureLocation(main_start_pos,main_end_pos,strand=strandToOutput)

			# 3. Create a SeqFeature
			main_feature = SeqFeature.SeqFeature(main_feature_location,type=main_feature_type, qualifiers=main_feature_qualifiers)
			'''
			#find d-loop part
			#basically just look for a big gap between last feature and this current feature, if there is a gap that
			#is about the size of a d-loop, it most likely is a dloop, since nothing aligned with it and it has that size
			#ignore this check if a cloroplast was built
			if lastFeatureAlignment != None and dLoopFound == False and buildCloroplast == False and dLoopSize > 0:
				if thisFeatureAlignment.startBase > lastFeatureAlignment.endBase + dLoopSize \
				 and thisFeatureAlignment.startBase < lastFeatureAlignment.endBase + 3200:
					dLoopFound = True
					dLoopStartPos = SeqFeature.ExactPosition(lastFeatureAlignment.endBase)
					dLoopEndPos = SeqFeature.ExactPosition(thisFeatureAlignment.startBase)
					dLoopLocation = SeqFeature.FeatureLocation(dLoopStartPos,dLoopEndPos,strand=-1)
					dLoopType = "D-loop"
					dLoopFeature = SeqFeature.SeqFeature(dLoopLocation,type=dLoopType)
					finalResults.features.append(dLoopFeature)

			lastFeatureAlignment = thisFeatureAlignment
			'''

			# 4. Append your newly created SeqFeature to your SeqRecord
			if main_feature_type == "gene":				
				cds_qualifiers = dict(main_feature_qualifiers)
				coding_dna = Seq(str(finalResults.seq[thisFeatureAlignment.startBase:thisFeatureAlignment.endBase]), IUPAC.unambiguous_dna)
				if strandToOutput == -1:
					coding_dna = coding_dna.reverse_complement()
				translationTable = thisFeatureAlignment.translationTable
				tableToUse = CodonTable.unambiguous_dna_by_id[translationTable]
				listOfStartCodons = []
				for startCodon in tableToUse.start_codons:
					startCodonSeq = Seq(startCodon, IUPAC.unambiguous_dna)
					startCodonTranslation = str(startCodonSeq.translate(table=translationTable))
					if startCodonTranslation not in listOfStartCodons:
						listOfStartCodons.append(startCodonTranslation)
					startCodons = tuple(listOfStartCodons) #need to make it a tuple so that startswith works with it
				stopCodons = ('*','$','#','+')
				nWalkStart = 20
				nWalkStop = 20
				'''
				For genes in the -1 strand, we look for the stop codons at the start and the start codons at the end!
				'''
				if strandToOutput == -1:
					tempStartCodons = startCodons
					tempStopCodons = stopCodons
					startCodons = tempStopCodons
					stopCodons = tempStartCodons
					nWalkStart = 20
					nWalkStop = 20			

				try:
					'''
					Making sure it starts with startCodons
					'''
					try:
						coding_dna_Translation = coding_dna.translate(table=translationTable)
						coding_dna_TranslationForward = coding_dna.translate(table=translationTable)
						coding_dna_TranslationBackward = coding_dna.translate(table=translationTable)
						if strandToOutput == -1:
							coding_dna_TranslationBackward = coding_dna_TranslationBackward[::-1]
							coding_dna_TranslationForward = coding_dna_TranslationForward[::-1]
						coding_dna_Forward = coding_dna
						coding_dna_Backward = coding_dna
						startBase = thisFeatureAlignment.startBase
						endBase = thisFeatureAlignment.endBase
						n = 0
						while not coding_dna_TranslationForward.startswith(startCodons) \
						and not coding_dna_TranslationBackward.startswith(startCodons) and n < nWalkStart and startBase - (3*n) >= 0:
							try:
								n += 1
								coding_dna_Backward = Seq(str(finalResults.seq[startBase - (3*n):endBase]), IUPAC.unambiguous_dna)
								if strandToOutput == -1:
									coding_dna_Backward = coding_dna_Backward.reverse_complement()
								coding_dna_TranslationBackward = coding_dna_Backward.translate(table=translationTable)
								coding_dna_Forward = Seq(str(finalResults.seq[startBase + (3*n):endBase]), IUPAC.unambiguous_dna)
								if strandToOutput == -1:
									coding_dna_Forward = coding_dna_Forward.reverse_complement()
								coding_dna_TranslationForward = coding_dna_Forward.translate(table=translationTable)
								if strandToOutput == -1:
									coding_dna_TranslationBackward = coding_dna_TranslationBackward[::-1]
									coding_dna_TranslationForward = coding_dna_TranslationForward[::-1]
							except:
								pass
						else:
							if coding_dna_TranslationBackward.startswith(startCodons):
								main_start_pos = SeqFeature.ExactPosition(startBase - (3*n))
								startBase += (3*n)
								thisFeatureAlignment.startBase = startBase
								main_feature_location = SeqFeature.FeatureLocation(main_start_pos,main_end_pos,strand=strandToOutput)
							elif coding_dna_TranslationForward.startswith(startCodons):
								main_start_pos = SeqFeature.ExactPosition(startBase + (3*n))
								startBase += (3*n)
								thisFeatureAlignment.startBase = startBase
								main_feature_location = SeqFeature.FeatureLocation(main_start_pos,main_end_pos,strand=strandToOutput)
					except:
						pass

					'''
					Making sure it ends with * (stop codon)
					'''
					try:
						coding_dna_Translation = coding_dna.translate(table=translationTable)
						coding_dna_TranslationForward = coding_dna.translate(table=translationTable)
						coding_dna_TranslationBackward = coding_dna.translate(table=translationTable)
						if strandToOutput == -1:
							coding_dna_TranslationBackward = coding_dna_TranslationBackward[::-1]
							coding_dna_TranslationForward = coding_dna_TranslationForward[::-1]
						coding_dna_Forward = coding_dna
						coding_dna_Backward = coding_dna
						startBase = thisFeatureAlignment.startBase
						endBase = thisFeatureAlignment.endBase
						n = 0
						while not coding_dna_TranslationForward.endswith(stopCodons) \
						and not coding_dna_TranslationBackward.endswith(stopCodons) and n < nWalkStop and endBase + (3*n) <= len(finalResults):
							try:
								n += 1
								coding_dna_Backward = Seq(str(finalResults.seq[startBase:endBase - (3*n)]), IUPAC.unambiguous_dna)
								if strandToOutput == -1:
									coding_dna_Backward = coding_dna_Backward.reverse_complement()
								coding_dna_TranslationBackward = coding_dna_Backward.translate(table=translationTable)
								
								coding_dna_Forward = Seq(str(finalResults.seq[startBase:endBase + (3*n)]), IUPAC.unambiguous_dna)
								if strandToOutput == -1:
									coding_dna_Forward = coding_dna_Forward.reverse_complement()
								coding_dna_TranslationForward = coding_dna_Forward.translate(table=translationTable)

								if strandToOutput == -1:
									coding_dna_TranslationBackward = coding_dna_TranslationBackward[::-1]
									coding_dna_TranslationForward = coding_dna_TranslationForward[::-1]
							except:
								pass
						else:
							if coding_dna_TranslationBackward.endswith(stopCodons):
								main_end_pos = SeqFeature.ExactPosition(endBase - (3 * n))
								endBase -= (3 * (n-1))
								thisFeatureAlignment.endBase = endBase
								main_feature_location = SeqFeature.FeatureLocation(main_start_pos,main_end_pos,strand=strandToOutput)
							elif coding_dna_TranslationForward.endswith(stopCodons):
								main_end_pos = SeqFeature.ExactPosition(endBase + (3 * n))
								endBase += (3 * (n-1))
								thisFeatureAlignment.endBase = endBase
								main_feature_location = SeqFeature.FeatureLocation(main_start_pos,main_end_pos,strand=strandToOutput)
					except:
						pass
					
					coding_dna = Seq(str(finalResults.seq[thisFeatureAlignment.startBase:thisFeatureAlignment.endBase]),IUPAC.unambiguous_dna)
					if strandToOutput == 1:
						coding_dna_Translation = coding_dna.translate(table=translationTable)
					else:
						coding_dna_Translation = coding_dna.reverse_complement().translate(table=translationTable)
					cds_qualifiers['translation'] = coding_dna_Translation
				except:
					cds_qualifiers['translation'] = 'ERROR'
				cds_feature = SeqFeature.SeqFeature(main_feature_location,type='CDS', qualifiers=cds_qualifiers)
				main_feature = SeqFeature.SeqFeature(main_feature_location,type=main_feature_type, qualifiers=main_feature_qualifiers)
				finalResults.features.append(main_feature)
				finalResults.features.append(cds_feature)
			else: #if it's a tRNA or rRNA
				gene_feature = SeqFeature.SeqFeature(main_feature_location,type='gene', qualifiers=main_feature_qualifiers)
				finalResults.features.append(gene_feature)
				finalResults.features.append(main_feature)

	#returns the final SeqRecord object, with all features, so that the script that called genbankOutput can output this result whatever way
	#it wants
	return finalResults
