#!/usr/bin/python
#Version: 1.2
#Authors: Allio Remi & Schomaker-Bastos Alex
#ISEM - CNRS - LAMPADA - IBQM - UFRJ

import sys
import os.path
import operator
import collections
from Bio import SeqIO, SeqFeature
from Bio.Alphabet import generic_dna, generic_protein, IUPAC
from Bio.Data import CodonTable
from Bio.Seq import Seq


def read_fasta(fp):
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield (name, ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name: yield (name, ''.join(seq))

file1=open(sys.argv[1])
fasta=open(sys.argv[1].split("_raw.gff")[0]+".fasta")
record = SeqIO.read(open(sys.argv[1].split("_raw.gff")[0]+".gb"), "genbank", generic_dna)

dico_start={}
dico_end={}

for feature in record.features:
	start=""
	stop=""
	if feature.type.lower() == 'cds':
		if 'gene' in feature.qualifiers:
			if 'translation' in feature.qualifiers:
					if feature.location.strand == -1:
						start=Seq(str(record.seq[feature.location.end-3 :feature.location.end]),IUPAC.unambiguous_dna)
						dico_start[feature.qualifiers['gene'][0]]=str(start.reverse_complement())
						stop=Seq(str(record.seq[feature.location.start :feature.location.start+3]),IUPAC.unambiguous_dna)
						dico_end[feature.qualifiers['gene'][0]]=str(stop.reverse_complement())
					else:
						dico_start[feature.qualifiers['gene'][0]]=str(Seq(str(record.seq[feature.location.start :feature.location.start+3]),IUPAC.unambiguous_dna))
						dico_end[feature.qualifiers['gene'][0]]=str(Seq(str(record.seq[feature.location.end-3 :feature.location.end]),IUPAC.unambiguous_dna))			
		if 'product' in feature.qualifiers:
			if 'translation' in feature.qualifiers:
					if feature.location.strand == -1:
						start=Seq(str(record.seq[feature.location.end-3 :feature.location.end]),IUPAC.unambiguous_dna)
						dico_start[feature.qualifiers['gene'][0]]=str(start.reverse_complement())
						stop=Seq(str(record.seq[feature.location.start :feature.location.start+3]),IUPAC.unambiguous_dna)
						dico_end[feature.qualifiers['gene'][0]]=str(stop.reverse_complement())
					else:
						dico_start[feature.qualifiers['gene'][0]]=str(Seq(str(record.seq[feature.location.start :feature.location.start+3]),IUPAC.unambiguous_dna))
						dico_end[feature.qualifiers['gene'][0]]=str(Seq(str(record.seq[feature.location.end-3 :feature.location.end]),IUPAC.unambiguous_dna))

tableToUse = CodonTable.unambiguous_dna_by_id[int(sys.argv[3])]
listOfStartCodons = []
listOfStopCodons = []

for startCodon in tableToUse.start_codons:
	if startCodon not in listOfStartCodons:
		listOfStartCodons.append(str(startCodon))
	startCodons = listOfStartCodons
	
for stopCodon in tableToUse.stop_codons:
	if stopCodon not in listOfStopCodons:
		listOfStopCodons.append(str(stopCodon))
	stopCodons = listOfStopCodons

seqID=sys.argv[2]

dico={}
dicotrna={}
dicof={}
c=0

for name, seq in read_fasta(fasta):
	length=len(seq)

for line in open(sys.argv[1]):
	dico[line.rstrip().split("\t")[8]]=line.rstrip()
	
for line in open(sys.argv[1]):
	t=1
	line=line.rstrip()
	chro=line.split("\t")[0]
	gene=line.split("\t")[8]
	start=int(line.split("\t")[3])
	stop=int(line.split("\t")[4])
	if "tRNA" in gene:
		for k, v in dico.items():
			v=v.rstrip()
			tchro=v.split("\t")[0]
			tgene=v.split("\t")[8]
			tstart=int(v.split("\t")[3])
			tstop=int(v.split("\t")[4])
			if tgene != gene and tchro == chro and tstart <= start and tstop >= start+6:
					t=0
					dicotrna[gene]=gene
					break
			elif tgene != gene and tchro == chro and tstart >= start and tstop <= stop:
					t=0
					dicotrna[gene]=gene
					break		
			elif tgene != gene and tchro == chro and tstart <= stop-6 and tstop >= stop:
					t=0
					dicotrna[gene]=gene
					break
	dicof[start]=line.rstrip()	
	
sorted_x = sorted(dicof.items(), key=operator.itemgetter(0))
sorted_dict = collections.OrderedDict(sorted_x)
gout=sys.argv[1].split("_raw.gff")[0]+".gff"
gout=open(gout, "w")
gout.write(seqID+"\t"+"mitofinder\tsource\t1\t"+str(length)+"\t.\t+\t.\tName=source "+seqID.replace("_"," ")+"\n")
tout=sys.argv[1].split("_raw.gff")[0]+".tbl"
tout=open(tout, "w")
tout.write(">Feature "+seqID+"\n")
#tout.write("1\t"+str(length)+"\tREFERENCE\n")
#tout.write("\t\t\tMitoFinder\txxxxxx\n")

for k, v in sorted_dict.items():
	if not dicotrna.has_key(v.split("\t")[8]):
		col1=seqID
		col2=v.split("\t")[1]
		col3=v.split("\t")[2]
		col4=v.split("\t")[3]
		col5=v.split("\t")[4]
		col6=v.split("\t")[5]
		col7=v.split("\t")[6]
		col8=v.split("\t")[7]
		col9=v.split("\t")[8].rstrip()
		ext=""
		if "COX1" == col9:
			size=int(col5)-int(col4)+1
			if size%3 == 0:
				if col7 == "+":
					if dico_start[col9] in startCodons:
						start=col4
					else:
						start="<"+col4
						ext="5' Partial CDS"
					if dico_end[col9] in stopCodons:
						stop=col5
					else:
						stop=">"+col5
						if ext != "":
							ext="Partial CDS"
						else:
							ext="3' Partial CDS"				
					tout.write(start+"\t"+stop+"\t"+"gene\n")
					tout.write("\t\t\tgene\t"+col9+"\n")
					tout.write(start+"\t"+stop+"\t"+"CDS\n")
					tout.write("\t\t\tproduct\tcytochrome c oxidase subunit I\n")
					tout.write("\t\t\ttransl_table\t"+sys.argv[3]+"\n")
					if ext != "":
						tout.write("\t\t\tnote\t"+ext+"\n")						
				if col7 == "-":
					if dico_start[col9] in startCodons:
						start=col5
					else:
						start="<"+col5
						ext="3' Partial CDS"
					if dico_end[col9] in stopCodons:
						stop=col4
					else:
						stop=">"+col4
						if ext != "":
							ext="Partial CDS"
						else:
							ext="5' Partial CDS"						
					tout.write(start+"\t"+stop+"\t"+"gene\n")
					tout.write("\t\t\tgene\t"+col9+"\n")
					tout.write(start+"\t"+stop+"\t"+"CDS\n")
					tout.write("\t\t\tproduct\tcytochrome c oxidase subunit I\n")
					tout.write("\t\t\ttransl_table\t"+sys.argv[3]+"\n")				
					if ext != "":
						tout.write("\t\t\tnote\t"+ext+"\n")						
			else:
				ext="Note: Not a multiple of 3"
				if col7 == "+":
					tout.write(col4+"\t"+">"+col5+"\t"+"gene\n")
					tout.write("\t\t\tgene\t"+col9+"\n")
					tout.write(col4+"\t"+">"+col5+"\t"+"CDS\n")
					tout.write("\t\t\tproduct\tcytochrome c oxidase subunit I\n")
					tout.write("\t\t\ttransl_table\t"+sys.argv[3]+"\n")
					tout.write("\t\t\tnote\tPartial sequence\n")
				if col7 == "-":
					tout.write("<"+col5+"\t"+col4+"\t"+"gene\n")
					tout.write("\t\t\tgene\t"+col9+"\n")
					tout.write("<"+col5+"\t"+col4+"\t"+"CDS\n")
					tout.write("\t\t\tproduct\tcytochrome c oxidase subunit I\n")
					tout.write("\t\t\ttransl_table\t"+sys.argv[3]+"\n")
					tout.write("\t\t\tnote\tPartial sequence\n")
			gout.write(col1+"\t"+col2+"\t"+"gene"+"\t"+col4+"\t"+col5+"\t"+col6+"\t"+col7+"\t"+col8+"\t"+"Name="+col9+" gene "+ext+"\n")
			gout.write(col1+"\t"+col2+"\t"+"CDS"+"\t"+col4+"\t"+col5+"\t"+col6+"\t"+col7+"\t"+col8+"\t"+"Name="+col9+" CDS "+ext+"\n")
		
		if "COX2" == col9:
			size=int(col5)-int(col4)+1
			if size%3 == 0:
				if col7 == "+":
					if dico_start[col9] in startCodons:
						start=col4
					else:
						start="<"+col4
						ext="5' Partial CDS"
					if dico_end[col9] in stopCodons:
						stop=col5
					else:
						stop=">"+col5
						if ext != "":
							ext="Partial CDS"
						else:
							ext="3' Partial CDS"				
					tout.write(start+"\t"+stop+"\t"+"gene\n")
					tout.write("\t\t\tgene\t"+col9+"\n")
					tout.write(start+"\t"+stop+"\t"+"CDS\n")
					tout.write("\t\t\tproduct\tcytochrome c oxidase subunit II\n")
					tout.write("\t\t\ttransl_table\t"+sys.argv[3]+"\n")
					if ext != "":
						tout.write("\t\t\tnote\t"+ext+"\n")						
				if col7 == "-":
					if dico_start[col9] in startCodons:
						start=col5
					else:
						start="<"+col5
						ext="3' Partial CDS"
					if dico_end[col9] in stopCodons:
						stop=col4
					else:
						stop=">"+col4
						if ext != "":
							ext="Partial CDS"
						else:
							ext="5' Partial CDS"						
					tout.write(start+"\t"+stop+"\t"+"gene\n")
					tout.write("\t\t\tgene\t"+col9+"\n")
					tout.write(start+"\t"+stop+"\t"+"CDS\n")
					tout.write("\t\t\tproduct\tcytochrome c oxidase subunit II\n")
					tout.write("\t\t\ttransl_table\t"+sys.argv[3]+"\n")				
					if ext != "":
						tout.write("\t\t\tnote\t"+ext+"\n")						
			else:
				ext="Note: Not a multiple of 3"
				if col7 == "+":
					tout.write(col4+"\t"+">"+col5+"\t"+"gene\n")
					tout.write("\t\t\tgene\t"+col9+"\n")
					tout.write(col4+"\t"+">"+col5+"\t"+"CDS\n")
					tout.write("\t\t\tproduct\tcytochrome c oxidase subunit II\n")
					tout.write("\t\t\ttransl_table\t"+sys.argv[3]+"\n")
					tout.write("\t\t\tnote\tPartial sequence\n")
				if col7 == "-":
					tout.write("<"+col5+"\t"+col4+"\t"+"gene\n")
					tout.write("\t\t\tgene\t"+col9+"\n")
					tout.write("<"+col5+"\t"+col4+"\t"+"CDS\n")
					tout.write("\t\t\tproduct\tcytochrome c oxidase subunit II\n")
					tout.write("\t\t\ttransl_table\t"+sys.argv[3]+"\n")
					tout.write("\t\t\tnote\tPartial sequence\n")
			gout.write(col1+"\t"+col2+"\t"+"gene"+"\t"+col4+"\t"+col5+"\t"+col6+"\t"+col7+"\t"+col8+"\t"+"Name="+col9+" gene "+ext+"\n")
			gout.write(col1+"\t"+col2+"\t"+"CDS"+"\t"+col4+"\t"+col5+"\t"+col6+"\t"+col7+"\t"+col8+"\t"+"Name="+col9+" CDS "+ext+"\n")
		
		if "COX3" == col9:
			size=int(col5)-int(col4)+1
			if size%3 == 0:
				if col7 == "+":
					if dico_start[col9] in startCodons:
						start=col4
					else:
						start="<"+col4
						ext="5' Partial CDS"
					if dico_end[col9] in stopCodons:
						stop=col5
					else:
						stop=">"+col5
						if ext != "":
							ext="Partial CDS"
						else:
							ext="3' Partial CDS"				
					tout.write(start+"\t"+stop+"\t"+"gene\n")
					tout.write("\t\t\tgene\t"+col9+"\n")
					tout.write(start+"\t"+stop+"\t"+"CDS\n")
					tout.write("\t\t\tproduct\tcytochrome c oxidase subunit III\n")
					tout.write("\t\t\ttransl_table\t"+sys.argv[3]+"\n")
					if ext != "":
						tout.write("\t\t\tnote\t"+ext+"\n")						
				if col7 == "-":
					if dico_start[col9] in startCodons:
						start=col5
					else:
						start="<"+col5
						ext="3' Partial CDS"
					if dico_end[col9] in stopCodons:
						stop=col4
					else:
						stop=">"+col4
						if ext != "":
							ext="Partial CDS"
						else:
							ext="5' Partial CDS"						
					tout.write(start+"\t"+stop+"\t"+"gene\n")
					tout.write("\t\t\tgene\t"+col9+"\n")
					tout.write(start+"\t"+stop+"\t"+"CDS\n")
					tout.write("\t\t\tproduct\tcytochrome c oxidase subunit III\n")
					tout.write("\t\t\ttransl_table\t"+sys.argv[3]+"\n")				
					if ext != "":
						tout.write("\t\t\tnote\t"+ext+"\n")						
			else:
				ext="Note: Not a multiple of 3"
				if col7 == "+":
					tout.write(col4+"\t"+">"+col5+"\t"+"gene\n")
					tout.write("\t\t\tgene\t"+col9+"\n")
					tout.write(col4+"\t"+">"+col5+"\t"+"CDS\n")
					tout.write("\t\t\tproduct\tcytochrome c oxidase subunit III\n")
					tout.write("\t\t\ttransl_table\t"+sys.argv[3]+"\n")
					tout.write("\t\t\tnote\tPartial sequence\n")
				if col7 == "-":
					tout.write("<"+col5+"\t"+col4+"\t"+"gene\n")
					tout.write("\t\t\tgene\t"+col9+"\n")
					tout.write("<"+col5+"\t"+col4+"\t"+"CDS\n")
					tout.write("\t\t\tproduct\tcytochrome c oxidase subunit III\n")
					tout.write("\t\t\ttransl_table\t"+sys.argv[3]+"\n")
					tout.write("\t\t\tnote\tPartial sequence\n")
			gout.write(col1+"\t"+col2+"\t"+"gene"+"\t"+col4+"\t"+col5+"\t"+col6+"\t"+col7+"\t"+col8+"\t"+"Name="+col9+" gene "+ext+"\n")
			gout.write(col1+"\t"+col2+"\t"+"CDS"+"\t"+col4+"\t"+col5+"\t"+col6+"\t"+col7+"\t"+col8+"\t"+"Name="+col9+" CDS "+ext+"\n")
	
		if "ND1" == col9:
			size=int(col5)-int(col4)+1
			if size%3 == 0:
				if col7 == "+":
					if dico_start[col9] in startCodons:
						start=col4
					else:
						start="<"+col4
						ext="5' Partial CDS"
					if dico_end[col9] in stopCodons:
						stop=col5
					else:
						stop=">"+col5
						if ext != "":
							ext="Partial CDS"
						else:
							ext="3' Partial CDS"				
					tout.write(start+"\t"+stop+"\t"+"gene\n")
					tout.write("\t\t\tgene\t"+col9+"\n")
					tout.write(start+"\t"+stop+"\t"+"CDS\n")
					tout.write("\t\t\tproduct\tNADH dehydrogenase subunit 1\n")
					tout.write("\t\t\ttransl_table\t"+sys.argv[3]+"\n")
					if ext != "":
						tout.write("\t\t\tnote\t"+ext+"\n")						
				if col7 == "-":
					if dico_start[col9] in startCodons:
						start=col5
					else:
						start="<"+col5
						ext="3' Partial CDS"
					if dico_end[col9] in stopCodons:
						stop=col4
					else:
						stop=">"+col4
						if ext != "":
							ext="Partial CDS"
						else:
							ext="5' Partial CDS"						
					tout.write(start+"\t"+stop+"\t"+"gene\n")
					tout.write("\t\t\tgene\t"+col9+"\n")
					tout.write(start+"\t"+stop+"\t"+"CDS\n")
					tout.write("\t\t\tproduct\tNADH dehydrogenase subunit 1\n")
					tout.write("\t\t\ttransl_table\t"+sys.argv[3]+"\n")				
					if ext != "":
						tout.write("\t\t\tnote\t"+ext+"\n")						
			else:
				ext="Note: Not a multiple of 3"
				if col7 == "+":
					tout.write(col4+"\t"+">"+col5+"\t"+"gene\n")
					tout.write("\t\t\tgene\t"+col9+"\n")
					tout.write(col4+"\t"+">"+col5+"\t"+"CDS\n")
					tout.write("\t\t\tproduct\tNADH dehydrogenase subunit 1\n")
					tout.write("\t\t\ttransl_table\t"+sys.argv[3]+"\n")
					tout.write("\t\t\tnote\tPartial sequence\n")
				if col7 == "-":
					tout.write("<"+col5+"\t"+col4+"\t"+"gene\n")
					tout.write("\t\t\tgene\t"+col9+"\n")
					tout.write("<"+col5+"\t"+col4+"\t"+"CDS\n")
					tout.write("\t\t\tproduct\tNADH dehydrogenase subunit 1\n")
					tout.write("\t\t\ttransl_table\t"+sys.argv[3]+"\n")
					tout.write("\t\t\tnote\tPartial sequence\n")
			gout.write(col1+"\t"+col2+"\t"+"gene"+"\t"+col4+"\t"+col5+"\t"+col6+"\t"+col7+"\t"+col8+"\t"+"Name="+col9+" gene "+ext+"\n")
			gout.write(col1+"\t"+col2+"\t"+"CDS"+"\t"+col4+"\t"+col5+"\t"+col6+"\t"+col7+"\t"+col8+"\t"+"Name="+col9+" CDS "+ext+"\n")

		if "ND2" == col9:
			size=int(col5)-int(col4)+1
			if size%3 == 0:
				if col7 == "+":
					if dico_start[col9] in startCodons:
						start=col4
					else:
						start="<"+col4
						ext="5' Partial CDS"
					if dico_end[col9] in stopCodons:
						stop=col5
					else:
						stop=">"+col5
						if ext != "":
							ext="Partial CDS"
						else:
							ext="3' Partial CDS"				
					tout.write(start+"\t"+stop+"\t"+"gene\n")
					tout.write("\t\t\tgene\t"+col9+"\n")
					tout.write(start+"\t"+stop+"\t"+"CDS\n")
					tout.write("\t\t\tproduct\tNADH dehydrogenase subunit 2\n")
					tout.write("\t\t\ttransl_table\t"+sys.argv[3]+"\n")
					if ext != "":
						tout.write("\t\t\tnote\t"+ext+"\n")						
				if col7 == "-":
					if dico_start[col9] in startCodons:
						start=col5
					else:
						start="<"+col5
						ext="3' Partial CDS"
					if dico_end[col9] in stopCodons:
						stop=col4
					else:
						stop=">"+col4
						if ext != "":
							ext="Partial CDS"
						else:
							ext="5' Partial CDS"						
					tout.write(start+"\t"+stop+"\t"+"gene\n")
					tout.write("\t\t\tgene\t"+col9+"\n")
					tout.write(start+"\t"+stop+"\t"+"CDS\n")
					tout.write("\t\t\tproduct\tNADH dehydrogenase subunit 2\n")
					tout.write("\t\t\ttransl_table\t"+sys.argv[3]+"\n")				
					if ext != "":
						tout.write("\t\t\tnote\t"+ext+"\n")						
			else:
				ext="Note: Not a multiple of 3"
				if col7 == "+":
					tout.write(col4+"\t"+">"+col5+"\t"+"gene\n")
					tout.write("\t\t\tgene\t"+col9+"\n")
					tout.write(col4+"\t"+">"+col5+"\t"+"CDS\n")
					tout.write("\t\t\tproduct\tNADH dehydrogenase subunit 2\n")
					tout.write("\t\t\ttransl_table\t"+sys.argv[3]+"\n")
					tout.write("\t\t\tnote\tPartial sequence\n")
				if col7 == "-":
					tout.write("<"+col5+"\t"+col4+"\t"+"gene\n")
					tout.write("\t\t\tgene\t"+col9+"\n")
					tout.write("<"+col5+"\t"+col4+"\t"+"CDS\n")
					tout.write("\t\t\tproduct\tNADH dehydrogenase subunit 2\n")
					tout.write("\t\t\ttransl_table\t"+sys.argv[3]+"\n")
					tout.write("\t\t\tnote\tPartial sequence\n")
			gout.write(col1+"\t"+col2+"\t"+"gene"+"\t"+col4+"\t"+col5+"\t"+col6+"\t"+col7+"\t"+col8+"\t"+"Name="+col9+" gene "+ext+"\n")
			gout.write(col1+"\t"+col2+"\t"+"CDS"+"\t"+col4+"\t"+col5+"\t"+col6+"\t"+col7+"\t"+col8+"\t"+"Name="+col9+" CDS "+ext+"\n")
		
		if "ND3" == col9:
			size=int(col5)-int(col4)+1
			if size%3 == 0:
				if col7 == "+":
					if dico_start[col9] in startCodons:
						start=col4
					else:
						start="<"+col4
						ext="5' Partial CDS"
					if dico_end[col9] in stopCodons:
						stop=col5
					else:
						stop=">"+col5
						if ext != "":
							ext="Partial CDS"
						else:
							ext="3' Partial CDS"				
					tout.write(start+"\t"+stop+"\t"+"gene\n")
					tout.write("\t\t\tgene\t"+col9+"\n")
					tout.write(start+"\t"+stop+"\t"+"CDS\n")
					tout.write("\t\t\tproduct\tNADH dehydrogenase subunit 3\n")
					tout.write("\t\t\ttransl_table\t"+sys.argv[3]+"\n")
					if ext != "":
						tout.write("\t\t\tnote\t"+ext+"\n")						
				if col7 == "-":
					if dico_start[col9] in startCodons:
						start=col5
					else:
						start="<"+col5
						ext="3' Partial CDS"
					if dico_end[col9] in stopCodons:
						stop=col4
					else:
						stop=">"+col4
						if ext != "":
							ext="Partial CDS"
						else:
							ext="5' Partial CDS"						
					tout.write(start+"\t"+stop+"\t"+"gene\n")
					tout.write("\t\t\tgene\t"+col9+"\n")
					tout.write(start+"\t"+stop+"\t"+"CDS\n")
					tout.write("\t\t\tproduct\tNADH dehydrogenase subunit 3\n")
					tout.write("\t\t\ttransl_table\t"+sys.argv[3]+"\n")				
					if ext != "":
						tout.write("\t\t\tnote\t"+ext+"\n")						
			else:
				ext="Note: Not a multiple of 3"
				if col7 == "+":
					tout.write(col4+"\t"+">"+col5+"\t"+"gene\n")
					tout.write("\t\t\tgene\t"+col9+"\n")
					tout.write(col4+"\t"+">"+col5+"\t"+"CDS\n")
					tout.write("\t\t\tproduct\tNADH dehydrogenase subunit 3\n")
					tout.write("\t\t\ttransl_table\t"+sys.argv[3]+"\n")
					tout.write("\t\t\tnote\tPartial sequence\n")
				if col7 == "-":
					tout.write("<"+col5+"\t"+col4+"\t"+"gene\n")
					tout.write("\t\t\tgene\t"+col9+"\n")
					tout.write("<"+col5+"\t"+col4+"\t"+"CDS\n")
					tout.write("\t\t\tproduct\tNADH dehydrogenase subunit 3\n")
					tout.write("\t\t\ttransl_table\t"+sys.argv[3]+"\n")
					tout.write("\t\t\tnote\tPartial sequence\n")
			gout.write(col1+"\t"+col2+"\t"+"gene"+"\t"+col4+"\t"+col5+"\t"+col6+"\t"+col7+"\t"+col8+"\t"+"Name="+col9+" gene "+ext+"\n")
			gout.write(col1+"\t"+col2+"\t"+"CDS"+"\t"+col4+"\t"+col5+"\t"+col6+"\t"+col7+"\t"+col8+"\t"+"Name="+col9+" CDS "+ext+"\n")

		if "ND4" == col9:
			size=int(col5)-int(col4)+1
			if size%3 == 0:
				if col7 == "+":
					if dico_start[col9] in startCodons:
						start=col4
					else:
						start="<"+col4
						ext="5' Partial CDS"
					if dico_end[col9] in stopCodons:
						stop=col5
					else:
						stop=">"+col5
						if ext != "":
							ext="Partial CDS"
						else:
							ext="3' Partial CDS"				
					tout.write(start+"\t"+stop+"\t"+"gene\n")
					tout.write("\t\t\tgene\t"+col9+"\n")
					tout.write(start+"\t"+stop+"\t"+"CDS\n")
					tout.write("\t\t\tproduct\tNADH dehydrogenase subunit 4\n")
					tout.write("\t\t\ttransl_table\t"+sys.argv[3]+"\n")
					if ext != "":
						tout.write("\t\t\tnote\t"+ext+"\n")						
				if col7 == "-":
					if dico_start[col9] in startCodons:
						start=col5
					else:
						start="<"+col5
						ext="3' Partial CDS"
					if dico_end[col9] in stopCodons:
						stop=col4
					else:
						stop=">"+col4
						if ext != "":
							ext="Partial CDS"
						else:
							ext="5' Partial CDS"						
					tout.write(start+"\t"+stop+"\t"+"gene\n")
					tout.write("\t\t\tgene\t"+col9+"\n")
					tout.write(start+"\t"+stop+"\t"+"CDS\n")
					tout.write("\t\t\tproduct\tNADH dehydrogenase subunit 4\n")
					tout.write("\t\t\ttransl_table\t"+sys.argv[3]+"\n")				
					if ext != "":
						tout.write("\t\t\tnote\t"+ext+"\n")						
			else:
				ext="Note: Not a multiple of 3"
				if col7 == "+":
					tout.write(col4+"\t"+">"+col5+"\t"+"gene\n")
					tout.write("\t\t\tgene\t"+col9+"\n")
					tout.write(col4+"\t"+">"+col5+"\t"+"CDS\n")
					tout.write("\t\t\tproduct\tNADH dehydrogenase subunit 4\n")
					tout.write("\t\t\ttransl_table\t"+sys.argv[3]+"\n")
					tout.write("\t\t\tnote\tPartial sequence\n")
				if col7 == "-":
					tout.write("<"+col5+"\t"+col4+"\t"+"gene\n")
					tout.write("\t\t\tgene\t"+col9+"\n")
					tout.write("<"+col5+"\t"+col4+"\t"+"CDS\n")
					tout.write("\t\t\tproduct\tNADH dehydrogenase subunit 4\n")
					tout.write("\t\t\ttransl_table\t"+sys.argv[3]+"\n")
					tout.write("\t\t\tnote\tPartial sequence\n")
			gout.write(col1+"\t"+col2+"\t"+"gene"+"\t"+col4+"\t"+col5+"\t"+col6+"\t"+col7+"\t"+col8+"\t"+"Name="+col9+" gene "+ext+"\n")
			gout.write(col1+"\t"+col2+"\t"+"CDS"+"\t"+col4+"\t"+col5+"\t"+col6+"\t"+col7+"\t"+col8+"\t"+"Name="+col9+" CDS "+ext+"\n")

		if "ND4L" == col9:
			size=int(col5)-int(col4)+1
			if size%3 == 0:
				if col7 == "+":
					if dico_start[col9] in startCodons:
						start=col4
					else:
						start="<"+col4
						ext="5' Partial CDS"
					if dico_end[col9] in stopCodons:
						stop=col5
					else:
						stop=">"+col5
						if ext != "":
							ext="Partial CDS"
						else:
							ext="3' Partial CDS"				
					tout.write(start+"\t"+stop+"\t"+"gene\n")
					tout.write("\t\t\tgene\t"+col9+"\n")
					tout.write(start+"\t"+stop+"\t"+"CDS\n")
					tout.write("\t\t\tproduct\tNADH dehydrogenase subunit 4L\n")
					tout.write("\t\t\ttransl_table\t"+sys.argv[3]+"\n")
					if ext != "":
						tout.write("\t\t\tnote\t"+ext+"\n")						
				if col7 == "-":
					if dico_start[col9] in startCodons:
						start=col5
					else:
						start="<"+col5
						ext="3' Partial CDS"
					if dico_end[col9] in stopCodons:
						stop=col4
					else:
						stop=">"+col4
						if ext != "":
							ext="Partial CDS"
						else:
							ext="5' Partial CDS"						
					tout.write(start+"\t"+stop+"\t"+"gene\n")
					tout.write("\t\t\tgene\t"+col9+"\n")
					tout.write(start+"\t"+stop+"\t"+"CDS\n")
					tout.write("\t\t\tproduct\tNADH dehydrogenase subunit 4L\n")
					tout.write("\t\t\ttransl_table\t"+sys.argv[3]+"\n")				
					if ext != "":
						tout.write("\t\t\tnote\t"+ext+"\n")						
			else:
				ext="Note: Not a multiple of 3"
				if col7 == "+":
					tout.write(col4+"\t"+">"+col5+"\t"+"gene\n")
					tout.write("\t\t\tgene\t"+col9+"\n")
					tout.write(col4+"\t"+">"+col5+"\t"+"CDS\n")
					tout.write("\t\t\tproduct\tNADH dehydrogenase subunit 4L\n")
					tout.write("\t\t\ttransl_table\t"+sys.argv[3]+"\n")
					tout.write("\t\t\tnote\tPartial sequence\n")
				if col7 == "-":
					tout.write("<"+col5+"\t"+col4+"\t"+"gene\n")
					tout.write("\t\t\tgene\t"+col9+"\n")
					tout.write("<"+col5+"\t"+col4+"\t"+"CDS\n")
					tout.write("\t\t\tproduct\tNADH dehydrogenase subunit 4L\n")
					tout.write("\t\t\ttransl_table\t"+sys.argv[3]+"\n")
					tout.write("\t\t\tnote\tPartial sequence\n")
			gout.write(col1+"\t"+col2+"\t"+"gene"+"\t"+col4+"\t"+col5+"\t"+col6+"\t"+col7+"\t"+col8+"\t"+"Name="+col9+" gene "+ext+"\n")
			gout.write(col1+"\t"+col2+"\t"+"CDS"+"\t"+col4+"\t"+col5+"\t"+col6+"\t"+col7+"\t"+col8+"\t"+"Name="+col9+" CDS "+ext+"\n")

		if "ND5" == col9:
			size=int(col5)-int(col4)+1
			if size%3 == 0:
				if col7 == "+":
					if dico_start[col9] in startCodons:
						start=col4
					else:
						start="<"+col4
						ext="5' Partial CDS"
					if dico_end[col9] in stopCodons:
						stop=col5
					else:
						stop=">"+col5
						if ext != "":
							ext="Partial CDS"
						else:
							ext="3' Partial CDS"				
					tout.write(start+"\t"+stop+"\t"+"gene\n")
					tout.write("\t\t\tgene\t"+col9+"\n")
					tout.write(start+"\t"+stop+"\t"+"CDS\n")
					tout.write("\t\t\tproduct\tNADH dehydrogenase subunit 5\n")
					tout.write("\t\t\ttransl_table\t"+sys.argv[3]+"\n")
					if ext != "":
						tout.write("\t\t\tnote\t"+ext+"\n")						
				if col7 == "-":
					if dico_start[col9] in startCodons:
						start=col5
					else:
						start="<"+col5
						ext="3' Partial CDS"
					if dico_end[col9] in stopCodons:
						stop=col4
					else:
						stop=">"+col4
						if ext != "":
							ext="Partial CDS"
						else:
							ext="5' Partial CDS"						
					tout.write(start+"\t"+stop+"\t"+"gene\n")
					tout.write("\t\t\tgene\t"+col9+"\n")
					tout.write(start+"\t"+stop+"\t"+"CDS\n")
					tout.write("\t\t\tproduct\tNADH dehydrogenase subunit 5\n")
					tout.write("\t\t\ttransl_table\t"+sys.argv[3]+"\n")				
					if ext != "":
						tout.write("\t\t\tnote\t"+ext+"\n")						
			else:
				ext="Note: Not a multiple of 3"
				if col7 == "+":
					tout.write(col4+"\t"+">"+col5+"\t"+"gene\n")
					tout.write("\t\t\tgene\t"+col9+"\n")
					tout.write(col4+"\t"+">"+col5+"\t"+"CDS\n")
					tout.write("\t\t\tproduct\tNADH dehydrogenase subunit 5\n")
					tout.write("\t\t\ttransl_table\t"+sys.argv[3]+"\n")
					tout.write("\t\t\tnote\tPartial sequence\n")
				if col7 == "-":
					tout.write("<"+col5+"\t"+col4+"\t"+"gene\n")
					tout.write("\t\t\tgene\t"+col9+"\n")
					tout.write("<"+col5+"\t"+col4+"\t"+"CDS\n")
					tout.write("\t\t\tproduct\tNADH dehydrogenase subunit 5\n")
					tout.write("\t\t\ttransl_table\t"+sys.argv[3]+"\n")
					tout.write("\t\t\tnote\tPartial sequence\n")
			gout.write(col1+"\t"+col2+"\t"+"gene"+"\t"+col4+"\t"+col5+"\t"+col6+"\t"+col7+"\t"+col8+"\t"+"Name="+col9+" gene "+ext+"\n")
			gout.write(col1+"\t"+col2+"\t"+"CDS"+"\t"+col4+"\t"+col5+"\t"+col6+"\t"+col7+"\t"+col8+"\t"+"Name="+col9+" CDS "+ext+"\n")

		if "CYTB" == col9:
			size=int(col5)-int(col4)+1
			if size%3 == 0:
				if col7 == "+":
					if dico_start[col9] in startCodons:
						start=col4
					else:
						start="<"+col4
						ext="5' Partial CDS"
					if dico_end[col9] in stopCodons:
						stop=col5
					else:
						stop=">"+col5
						if ext != "":
							ext="Partial CDS"
						else:
							ext="3' Partial CDS"				
					tout.write(start+"\t"+stop+"\t"+"gene\n")
					tout.write("\t\t\tgene\t"+col9+"\n")
					tout.write(start+"\t"+stop+"\t"+"CDS\n")
					tout.write("\t\t\tproduct\tcytochrome b\n")
					tout.write("\t\t\ttransl_table\t"+sys.argv[3]+"\n")
					if ext != "":
						tout.write("\t\t\tnote\t"+ext+"\n")						
				if col7 == "-":
					if dico_start[col9] in startCodons:
						start=col5
					else:
						start="<"+col5
						ext="3' Partial CDS"
					if dico_end[col9] in stopCodons:
						stop=col4
					else:
						stop=">"+col4
						if ext != "":
							ext="Partial CDS"
						else:
							ext="5' Partial CDS"						
					tout.write(start+"\t"+stop+"\t"+"gene\n")
					tout.write("\t\t\tgene\t"+col9+"\n")
					tout.write(start+"\t"+stop+"\t"+"CDS\n")
					tout.write("\t\t\tproduct\tcytochrome b\n")
					tout.write("\t\t\ttransl_table\t"+sys.argv[3]+"\n")				
					if ext != "":
						tout.write("\t\t\tnote\t"+ext+"\n")						
			else:
				ext="Note: Not a multiple of 3"
				if col7 == "+":
					tout.write(col4+"\t"+">"+col5+"\t"+"gene\n")
					tout.write("\t\t\tgene\t"+col9+"\n")
					tout.write(col4+"\t"+">"+col5+"\t"+"CDS\n")
					tout.write("\t\t\tproduct\tcytochrome b\n")
					tout.write("\t\t\ttransl_table\t"+sys.argv[3]+"\n")
					tout.write("\t\t\tnote\tPartial sequence\n")
				if col7 == "-":
					tout.write("<"+col5+"\t"+col4+"\t"+"gene\n")
					tout.write("\t\t\tgene\t"+col9+"\n")
					tout.write("<"+col5+"\t"+col4+"\t"+"CDS\n")
					tout.write("\t\t\tproduct\tcytochrome b\n")
					tout.write("\t\t\ttransl_table\t"+sys.argv[3]+"\n")
					tout.write("\t\t\tnote\tPartial sequence\n")
			gout.write(col1+"\t"+col2+"\t"+"gene"+"\t"+col4+"\t"+col5+"\t"+col6+"\t"+col7+"\t"+col8+"\t"+"Name="+col9+" gene "+ext+"\n")
			gout.write(col1+"\t"+col2+"\t"+"CDS"+"\t"+col4+"\t"+col5+"\t"+col6+"\t"+col7+"\t"+col8+"\t"+"Name="+col9+" CDS "+ext+"\n")

		if "ATP8" == col9:
			size=int(col5)-int(col4)+1
			if size%3 == 0:
				if col7 == "+":
					if dico_start[col9] in startCodons:
						start=col4
					else:
						start="<"+col4
						ext="5' Partial CDS"
					if dico_end[col9] in stopCodons:
						stop=col5
					else:
						stop=">"+col5
						if ext != "":
							ext="Partial CDS"
						else:
							ext="3' Partial CDS"				
					tout.write(start+"\t"+stop+"\t"+"gene\n")
					tout.write("\t\t\tgene\t"+col9+"\n")
					tout.write(start+"\t"+stop+"\t"+"CDS\n")
					tout.write("\t\t\tproduct\tATP synthase F0 subunit 8\n")
					tout.write("\t\t\ttransl_table\t"+sys.argv[3]+"\n")
					if ext != "":
						tout.write("\t\t\tnote\t"+ext+"\n")						
				if col7 == "-":
					if dico_start[col9] in startCodons:
						start=col5
					else:
						start="<"+col5
						ext="3' Partial CDS"
					if dico_end[col9] in stopCodons:
						stop=col4
					else:
						stop=">"+col4
						if ext != "":
							ext="Partial CDS"
						else:
							ext="5' Partial CDS"						
					tout.write(start+"\t"+stop+"\t"+"gene\n")
					tout.write("\t\t\tgene\t"+col9+"\n")
					tout.write(start+"\t"+stop+"\t"+"CDS\n")
					tout.write("\t\t\tproduct\tATP synthase F0 subunit 8\n")
					tout.write("\t\t\ttransl_table\t"+sys.argv[3]+"\n")	
					if ext != "":
						tout.write("\t\t\tnote\t"+ext+"\n")						
			else:
				ext="Note: Not a multiple of 3"
				if col7 == "+":
					tout.write(col4+"\t"+">"+col5+"\t"+"gene\n")
					tout.write("\t\t\tgene\t"+col9+"\n")
					tout.write(col4+"\t"+">"+col5+"\t"+"CDS\n")
					tout.write("\t\t\tproduct\tATP synthase F0 subunit 8\n")
					tout.write("\t\t\ttransl_table\t"+sys.argv[3]+"\n")
					tout.write("\t\t\tnote\tPartial sequence\n")
				if col7 == "-":
					tout.write("<"+col5+"\t"+col4+"\t"+"gene\n")
					tout.write("\t\t\tgene\t"+col9+"\n")
					tout.write("<"+col5+"\t"+col4+"\t"+"CDS\n")
					tout.write("\t\t\tproduct\tATP synthase F0 subunit 8\n")
					tout.write("\t\t\ttransl_table\t"+sys.argv[3]+"\n")
					tout.write("\t\t\tnote\tPartial sequence\n")
			gout.write(col1+"\t"+col2+"\t"+"gene"+"\t"+col4+"\t"+col5+"\t"+col6+"\t"+col7+"\t"+col8+"\t"+"Name="+col9+" gene "+ext+"\n")
			gout.write(col1+"\t"+col2+"\t"+"CDS"+"\t"+col4+"\t"+col5+"\t"+col6+"\t"+col7+"\t"+col8+"\t"+"Name="+col9+" CDS "+ext+"\n")

		if "ATP6" == col9:
			size=int(col5)-int(col4)+1
			if size%3 == 0:
				if col7 == "+":
					if dico_start[col9] in startCodons:
						start=col4
					else:
						start="<"+col4
						ext="5' Partial CDS"
					if dico_end[col9] in stopCodons:
						stop=col5
					else:
						stop=">"+col5
						if ext != "":
							ext="Partial CDS"
						else:
							ext="3' Partial CDS"				
					tout.write(start+"\t"+stop+"\t"+"gene\n")
					tout.write("\t\t\tgene\t"+col9+"\n")
					tout.write(start+"\t"+stop+"\t"+"CDS\n")
					tout.write("\t\t\tproduct\tATP synthase F0 subunit 6\n")
					tout.write("\t\t\ttransl_table\t"+sys.argv[3]+"\n")
					if ext != "":
						tout.write("\t\t\tnote\t"+ext+"\n")						
				if col7 == "-":
					if dico_start[col9] in startCodons:
						start=col5
					else:
						start="<"+col5
						ext="3' Partial CDS"
					if dico_end[col9] in stopCodons:
						stop=col4
					else:
						stop=">"+col4
						if ext != "":
							ext="Partial CDS"
						else:
							ext="5' Partial CDS"						
					tout.write(start+"\t"+stop+"\t"+"gene\n")
					tout.write("\t\t\tgene\t"+col9+"\n")
					tout.write(start+"\t"+stop+"\t"+"CDS\n")
					tout.write("\t\t\tproduct\tATP synthase F0 subunit 6\n")
					tout.write("\t\t\ttransl_table\t"+sys.argv[3]+"\n")	
					if ext != "":
						tout.write("\t\t\tnote\t"+ext+"\n")						
			else:
				ext="Note: Not a multiple of 3"
				if col7 == "+":
					tout.write(col4+"\t"+">"+col5+"\t"+"gene\n")
					tout.write("\t\t\tgene\t"+col9+"\n")
					tout.write(col4+"\t"+">"+col5+"\t"+"CDS\n")
					tout.write("\t\t\tproduct\tATP synthase F0 subunit 6\n")
					tout.write("\t\t\ttransl_table\t"+sys.argv[3]+"\n")
					tout.write("\t\t\tnote\tPartial sequence\n")
				if col7 == "-":
					tout.write("<"+col5+"\t"+col4+"\t"+"gene\n")
					tout.write("\t\t\tgene\t"+col9+"\n")
					tout.write("<"+col5+"\t"+col4+"\t"+"CDS\n")
					tout.write("\t\t\tproduct\tATP synthase F0 subunit 6\n")
					tout.write("\t\t\ttransl_table\t"+sys.argv[3]+"\n")
					tout.write("\t\t\tnote\tPartial sequence\n")
			gout.write(col1+"\t"+col2+"\t"+"gene"+"\t"+col4+"\t"+col5+"\t"+col6+"\t"+col7+"\t"+col8+"\t"+"Name="+col9+" gene "+ext+"\n")
			gout.write(col1+"\t"+col2+"\t"+"CDS"+"\t"+col4+"\t"+col5+"\t"+col6+"\t"+col7+"\t"+col8+"\t"+"Name="+col9+" CDS "+ext+"\n")

		if "tRNA" in col9:
			size=int(col5)-int(col4)+1
			if size <= 120:
				if col9 == "tRNA-???":
					col9 = "tRNA-Xxx"
				if col7 == "+":
					tout.write(col4+"\t"+col5+"\t"+"gene\n")
					tout.write("\t\t\tgene\t"+col9+"\n")
					tout.write(col4+"\t"+col5+"\t"+"tRNA\n")
					tout.write("\t\t\tproduct\t"+col9+"\n")
				if col7 == "-":
					tout.write(col5+"\t"+col4+"\t"+"gene\n")
					tout.write("\t\t\tgene\t"+col9+"\n")
					tout.write(col5+"\t"+col4+"\t"+"tRNA\n")
					tout.write("\t\t\tproduct\t"+col9+"\n")					
				gout.write(col1+"\t"+col2+"\t"+"gene"+"\t"+col4+"\t"+col5+"\t"+col6+"\t"+col7+"\t"+col8+"\t"+"Name="+col9+" gene"+"\n")
				gout.write(col1+"\t"+col2+"\t"+"tRNA"+"\t"+col4+"\t"+col5+"\t"+col6+"\t"+col7+"\t"+col8+"\t"+"Name="+col9+" "+"\n")
		if "trn" in col9:
			size=int(col5)-int(col4)+1
			if size <= 120:
				if col9 == "tRNA-???":
					col9 = "tRNA-Xxx"
				if col7 == "+":
					tout.write(col4+"\t"+col5+"\t"+"gene\n")
					tout.write("\t\t\tgene\t"+col9+"\n")
					tout.write(col4+"\t"+col5+"\t"+"tRNA\n")
					tout.write("\t\t\tproduct\t"+col9+"\n")
				if col7 == "-":
					tout.write(col5+"\t"+col4+"\t"+"gene\n")
					tout.write("\t\t\tgene\t"+col9+"\n")
					tout.write(col5+"\t"+col4+"\t"+"tRNA\n")
					tout.write("\t\t\tproduct\t"+col9+"\n")	
				gout.write(col1+"\t"+col2+"\t"+"gene"+"\t"+col4+"\t"+col5+"\t"+col6+"\t"+col7+"\t"+col8+"\t"+"Name="+col9+" gene"+"\n")
				gout.write(col1+"\t"+col2+"\t"+"tRNA"+"\t"+col4+"\t"+col5+"\t"+col6+"\t"+col7+"\t"+col8+"\t"+"Name="+col9+" tRNA"+"\n")
		if "rrnL" in col9:
			if col7 == "+":
				tout.write(col4+"\t"+col5+"\t"+"gene\n")
				tout.write("\t\t\tgene\t"+col9+"\n")
				tout.write(col4+"\t"+col5+"\t"+"rRNA\n")
				tout.write("\t\t\tgene\t"+col9+"\n")
				tout.write("\t\t\tproduct\t16S ribosomal RNA\n")
			if col7 == "-":
				tout.write(col5+"\t"+col4+"\t"+"gene\n")
				tout.write("\t\t\tgene\t"+col9+"\n")
				tout.write(col5+"\t"+col4+"\t"+"rRNA\n")
				tout.write("\t\t\tgene\t"+col9+"\n")
				tout.write("\t\t\tproduct\t16S ribosomal RNA\n")				
			gout.write(col1+"\t"+col2+"\t"+"gene"+"\t"+col4+"\t"+col5+"\t"+col6+"\t"+col7+"\t"+col8+"\t"+"Name="+col9+" gene"+"\n")
			gout.write(col1+"\t"+col2+"\t"+"rRNA"+"\t"+col4+"\t"+col5+"\t"+col6+"\t"+col7+"\t"+col8+"\t"+"Name="+col9+" rRNA"+"\n")
		if "rrnS" in col9:
			if col7 == "+":
				tout.write(col4+"\t"+col5+"\t"+"gene\n")
				tout.write("\t\t\tgene\t"+col9+"\n")
				tout.write(col4+"\t"+col5+"\t"+"rRNA\n")
				tout.write("\t\t\tgene\t"+col9+"\n")
				tout.write("\t\t\tproduct\t12S ribosomal RNA\n")
			if col7 == "-":
				tout.write(col5+"\t"+col4+"\t"+"gene\n")
				tout.write("\t\t\tgene\t"+col9+"\n")
				tout.write(col5+"\t"+col4+"\t"+"rRNA\n")
				tout.write("\t\t\tgene\t"+col9+"\n")
				tout.write("\t\t\tproduct\t12S ribosomal RNA\n")
			gout.write(col1+"\t"+col2+"\t"+"gene"+"\t"+col4+"\t"+col5+"\t"+col6+"\t"+col7+"\t"+col8+"\t"+"Name="+col9+" gene"+"\n")
			gout.write(col1+"\t"+col2+"\t"+"rRNA"+"\t"+col4+"\t"+col5+"\t"+col6+"\t"+col7+"\t"+col8+"\t"+"Name="+col9+" rRNA"+"\n")
gout.close()
tout.close()

	
	
	
