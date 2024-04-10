#!/usr/bin/python2.7
#Version: 1.4
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
dico_lstart={}
dico_lend={}

for feature in record.features:
	start=""
	stop=""
	if feature.type.lower() == 'cds':
		if 'gene' in feature.qualifiers:
			if 'translation' in feature.qualifiers:
					if feature.location.strand == -1:
						start=Seq(str(record.seq[feature.location.end-3 :feature.location.end]),IUPAC.unambiguous_dna)
						
						if not dico_lstart.has_key(feature.qualifiers['gene'][0].split("_")[0]):
							dico_lstart[feature.qualifiers['gene'][0].split("_")[0]]=feature.location.end
							dico_start[feature.qualifiers['gene'][0].split("_")[0]]=str(start.reverse_complement())
						else:
							if dico_lstart.get(feature.qualifiers['gene'][0].split("_")[0]) < feature.location.end:
								dico_lstart[feature.qualifiers['gene'][0].split("_")[0]]=feature.location.end
								dico_start[feature.qualifiers['gene'][0].split("_")[0]]=str(start.reverse_complement())
						stop=Seq(str(record.seq[feature.location.start :feature.location.start+3]),IUPAC.unambiguous_dna)
						if not dico_lend.has_key(feature.qualifiers['gene'][0].split("_")[0]):
							dico_lend[feature.qualifiers['gene'][0].split("_")[0]]=feature.location.start
							dico_end[feature.qualifiers['gene'][0].split("_")[0]]=str(stop.reverse_complement())
						else:
							if dico_lend.get(feature.qualifiers['gene'][0].split("_")[0]) > feature.location.start:
								dico_lend[feature.qualifiers['gene'][0].split("_")[0]]=feature.location.start
								dico_end[feature.qualifiers['gene'][0].split("_")[0]]=str(stop.reverse_complement())
					else:
						if not dico_lstart.has_key(feature.qualifiers['gene'][0].split("_")[0]):
							dico_lstart[feature.qualifiers['gene'][0].split("_")[0]]=feature.location.start
							dico_start[feature.qualifiers['gene'][0].split("_")[0]]=str(Seq(str(record.seq[feature.location.start :feature.location.start+3]),IUPAC.unambiguous_dna))
						else:
							if dico_lstart.get(feature.qualifiers['gene'][0].split("_")[0]) > feature.location.start:
								dico_lstart[feature.qualifiers['gene'][0].split("_")[0]]=feature.location.start
								dico_start[feature.qualifiers['gene'][0].split("_")[0]]=str(Seq(str(record.seq[feature.location.start :feature.location.start+3]),IUPAC.unambiguous_dna))
						if not dico_lend.has_key(feature.qualifiers['gene'][0].split("_")[0]):
							dico_lend[feature.qualifiers['gene'][0].split("_")[0]]=feature.location.end
							dico_end[feature.qualifiers['gene'][0].split("_")[0]]=str(Seq(str(record.seq[feature.location.end-3 :feature.location.end]),IUPAC.unambiguous_dna))
						else:
							if dico_lend.get(feature.qualifiers['gene'][0].split("_")[0]) < feature.location.end:
								dico_lend[feature.qualifiers['gene'][0].split("_")[0]]=feature.location.end
								dico_end[feature.qualifiers['gene'][0].split("_")[0]]=str(Seq(str(record.seq[feature.location.end-3 :feature.location.end]),IUPAC.unambiguous_dna))
		if 'product' in feature.qualifiers:
					if feature.location.strand == -1:
						start=Seq(str(record.seq[feature.location.end-3 :feature.location.end]),IUPAC.unambiguous_dna)
						
						if not dico_lstart.has_key(feature.qualifiers['gene'][0].split("_")[0]):
							dico_lstart[feature.qualifiers['gene'][0].split("_")[0]]=feature.location.end
							dico_start[feature.qualifiers['gene'][0].split("_")[0]]=str(start.reverse_complement())
						else:
							if dico_lstart.get(feature.qualifiers['gene'][0].split("_")[0]) < feature.location.end:
								dico_lstart[feature.qualifiers['gene'][0].split("_")[0]]=feature.location.end
								dico_start[feature.qualifiers['gene'][0].split("_")[0]]=str(start.reverse_complement())
						stop=Seq(str(record.seq[feature.location.start :feature.location.start+3]),IUPAC.unambiguous_dna)
						if not dico_lend.has_key(feature.qualifiers['gene'][0].split("_")[0]):
							dico_lend[feature.qualifiers['gene'][0].split("_")[0]]=feature.location.start
							dico_end[feature.qualifiers['gene'][0].split("_")[0]]=str(stop.reverse_complement())
						else:
							if dico_lend.get(feature.qualifiers['gene'][0].split("_")[0]) > feature.location.start:
								dico_lend[feature.qualifiers['gene'][0].split("_")[0]]=feature.location.start
								dico_end[feature.qualifiers['gene'][0].split("_")[0]]=str(stop.reverse_complement())
					else:
						if not dico_lstart.has_key(feature.qualifiers['gene'][0].split("_")[0]):
							dico_lstart[feature.qualifiers['gene'][0].split("_")[0]]=feature.location.start
							dico_start[feature.qualifiers['gene'][0].split("_")[0]]=str(Seq(str(record.seq[feature.location.start :feature.location.start+3]),IUPAC.unambiguous_dna))
						else:
							if dico_lstart.get(feature.qualifiers['gene'][0].split("_")[0]) > feature.location.start:
								dico_lstart[feature.qualifiers['gene'][0].split("_")[0]]=feature.location.start
								dico_start[feature.qualifiers['gene'][0].split("_")[0]]=str(Seq(str(record.seq[feature.location.start :feature.location.start+3]),IUPAC.unambiguous_dna))
						if not dico_lend.has_key(feature.qualifiers['gene'][0].split("_")[0]):
							dico_lend[feature.qualifiers['gene'][0].split("_")[0]]=feature.location.end
							dico_end[feature.qualifiers['gene'][0].split("_")[0]]=str(Seq(str(record.seq[feature.location.end-3 :feature.location.end]),IUPAC.unambiguous_dna))
						else:
							if dico_lend.get(feature.qualifiers['gene'][0].split("_")[0]) < feature.location.end:
								dico_lend[feature.qualifiers['gene'][0].split("_")[0]]=feature.location.end
								dico_end[feature.qualifiers['gene'][0].split("_")[0]]=str(Seq(str(record.seq[feature.location.end-3 :feature.location.end]),IUPAC.unambiguous_dna))
		"""print feature.qualifiers['gene'][0].split("_")[0]
		print dico_end[feature.qualifiers['gene'][0].split("_")[0]]
		print dico_lend[feature.qualifiers['gene'][0].split("_")[0]]"""

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
"""print stopCodons
print startCodons"""
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
			if tgene != gene and tchro == chro and tstart <= start and tstop >= stop:
					print gene+" = "+str(start)+","+str(stop)
					print tgene+" = "+str(tstart)+","+str(tstop)
					t=0
					dicotrna[gene]=gene
					break
			elif tgene != gene and tchro == chro and tstart >= start and tstop <= stop:
					print gene+" = "+str(start)+","+str(stop)
					print tgene+" = "+str(tstart)+","+str(tstop)
					t=0
					dicotrna[gene]=gene
					break		
			"""elif tgene != gene and tchro == chro and tstart <= start and tstop >= stop:
					t=0
					dicotrna[gene]=gene
					break"""
	dicof[start]=line.rstrip()	

dicog={}
dicogl={}
for line in open(sys.argv[1]):
	line=line.rstrip()
	gene=line.split("\t")[8]
	if not "_" in gene:
		dicog[gene]=1
		dicogl[gene]=line
		
for line in open(sys.argv[1]):
	line=line.rstrip()
	gene=line.split("\t")[8]
	if sys.argv[4] == "no":
		seqID=line.split("\t")[0]
	if "_" in gene:
		if line.split("\t")[6]== "-":
			if int(line.split("\t")[3]) > int(dicogl.get(gene.split("_")[0]).split(";")[0].split("\t")[3]): 
				dicog[gene.split("_")[0]]+=1
				dicogl[gene.split("_")[0]]=line+";"+dicogl.get(gene.split("_")[0])
			else:
				dicog[gene.split("_")[0]]+=1
				dicogl[gene.split("_")[0]]=dicogl.get(gene.split("_")[0])+";"+line
		elif line.split("\t")[6]== "+":
			if int(line.split("\t")[3]) < int(dicogl.get(gene.split("_")[0]).split(";")[0].split("\t")[3]): 
				dicog[gene.split("_")[0]]+=1
				dicogl[gene.split("_")[0]]=line+";"+dicogl.get(gene.split("_")[0])
			else:
				dicog[gene.split("_")[0]]+=1
				dicogl[gene.split("_")[0]]=dicogl.get(gene.split("_")[0])+";"+line


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

dico_product={}
dico_product["COX1"]="cytochrome c oxidase subunit I"
dico_product["COX2"]="cytochrome c oxidase subunit II"
dico_product["COX3"]="cytochrome c oxidase subunit III"
dico_product["ND1"]="NADH dehydrogenase subunit 1"
dico_product["ND2"]="NADH dehydrogenase subunit 2"
dico_product["ND3"]="NADH dehydrogenase subunit 3"
dico_product["ND4"]="NADH dehydrogenase subunit 4"
dico_product["ND4L"]="NADH dehydrogenase subunit 4L"
dico_product["ND5"]="NADH dehydrogenase subunit 5"
dico_product["ND6"]="NADH dehydrogenase subunit 6"
dico_product["CYTB"]="cytochrome b"
dico_product["ATP6"]="ATP synthase F0 subunit 6"
dico_product["ATP8"]="ATP synthase F0 subunit 8"
dico_product["rrnL"]="16S ribosomal RNA"
dico_product["rrnS"]="12S ribosomal RNA"

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
		if not "trn" in col9 and not "tRNA" in col9 and not "rrn" in col9:
			if dicog.get(col9) == 1:
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
						if dico_product.has_key(col9):
							tout.write("\t\t\tproduct\t"+dico_product.get(col9)+"\n")
						else:
							tout.write("\t\t\tproduct\tunknown\n")
						tout.write("\t\t\ttransl_table\t"+sys.argv[3]+"\n")
						if ext != "":
							tout.write("\t\t\tnote\t"+ext+"\n")						
					if col7 == "-":
						if dico_start[col9] in startCodons:
							start=col5
						else:
							start="<"+col5
							ext="5' Partial CDS"
						if dico_end[col9] in stopCodons:
							stop=col4
						else:
							stop=">"+col4
							if ext != "":
								ext="Partial CDS"
							else:
								ext="3' Partial CDS"						
						tout.write(start+"\t"+stop+"\t"+"gene\n")
						tout.write("\t\t\tgene\t"+col9+"\n")
						tout.write(start+"\t"+stop+"\t"+"CDS\n")
						if dico_product.has_key(col9):
							tout.write("\t\t\tproduct\t"+dico_product.get(col9)+"\n")
						else:
							tout.write("\t\t\tproduct\t"+col9+"\n")
						tout.write("\t\t\ttransl_table\t"+sys.argv[3]+"\n")				
						if ext != "":
							tout.write("\t\t\tnote\t"+ext+"\n")						
				else:
					ext="Note: Not a multiple of 3"
					if col7 == "+":
						tout.write(col4+"\t"+">"+col5+"\t"+"gene\n")
						tout.write("\t\t\tgene\t"+col9+"\n")
						tout.write(col4+"\t"+">"+col5+"\t"+"CDS\n")
						tout.write("\t\t\tproduct\t"+dico_product.get(col9)+"\n")
						tout.write("\t\t\ttransl_table\t"+sys.argv[3]+"\n")
						tout.write("\t\t\tnote\tPartial sequence\n")
					if col7 == "-":
						tout.write("<"+col5+"\t"+col4+"\t"+"gene\n")
						tout.write("\t\t\tgene\t"+col9+"\n")
						tout.write("<"+col5+"\t"+col4+"\t"+"CDS\n")
						tout.write("\t\t\tproduct\t"+dico_product.get(col9)+"\n")
						tout.write("\t\t\ttransl_table\t"+sys.argv[3]+"\n")
						tout.write("\t\t\tnote\tPartial sequence\n")
				gout.write(col1+"\t"+col2+"\t"+"gene"+"\t"+col4+"\t"+col5+"\t"+col6+"\t"+col7+"\t"+col8+"\t"+"Name="+col9+" gene "+ext+"\n")
				gout.write(col1+"\t"+col2+"\t"+"CDS"+"\t"+col4+"\t"+col5+"\t"+col6+"\t"+col7+"\t"+col8+"\t"+"Name="+col9+" CDS "+ext+"\n")
			
			if dicog.get(col9) > 1:
				c=1
				list_tmp_gene=[]
				list_tmp_cds=[]
				dico_tmp_gff={}
				for line in dicogl.get(col9).split(";"):
					col1=seqID
					col2=line.split("\t")[1]
					col3=line.split("\t")[2]
					col4=line.split("\t")[3]
					col5=line.split("\t")[4]
					col6=line.split("\t")[5]
					col7=line.split("\t")[6]
					col8=line.split("\t")[7]
					col9=line.split("\t")[8].split("_")[0].rstrip()
					size=int(col5)-int(col4)+1
					if c == 1:
						if size%3 == 0:
							if col7 == "+":
								if dico_start[col9] in startCodons:
									start=col4
								else:
									start="<"+col4
									ext="5' Partial CDS"
								stop=col5
							if col7 == "-":
								if dico_end[col9] in stopCodons:
									stop=col4
								else:
									stop=">"+col4
									if ext != "":
										ext="Partial CDS"
									else:
										ext="3' Partial CDS"						
								start=col5
							list_tmp_gene.append(int(start.replace("<","").replace(">","")))
							list_tmp_gene.append(int(stop.replace("<","").replace(">","")))
							list_tmp_cds.append(start+"\t"+stop+"\t"+"CDS")
							dico_tmp_gff[int(col4)+1]=col1+"\t"+col2+"\t"+"CDS"+"\t"+col4+"\t"+col5+"\t"+col6+"\t"+col7+"\t"+col8+"\t"+"Parent="+col9+".1 CDS "+ext
							
					if c > 1 and c < int(dicog.get(col9)):
						if size%3 == 0:
							if col7 == "+":
								start=col4
								stop=col5
							if col7 == "-":
								start=col5
								stop=col4
							col9=line.split("\t")[8].split("_")[0].rstrip()
							list_tmp_gene.append(int(start.replace("<","").replace(">","")))
							list_tmp_gene.append(int(stop.replace("<","").replace(">","")))
							list_tmp_cds.append(start+"\t"+stop)
							dico_tmp_gff[int(col4)+1]=col1+"\t"+col2+"\t"+"CDS"+"\t"+col4+"\t"+col5+"\t"+col6+"\t"+col7+"\t"+col8+"\t"+"Parent="+col9+".1 CDS "+ext
							
					if c > 1 and c == int(dicog.get(col9)):
						if size%3 == 0:
							if col7 == "+":
								start=col4
								if dico_end[col9] in stopCodons:
									stop=col5
								else:
									stop=">"+col5
									if ext != "":
										ext="Partial CDS"
									else:
										ext="3' Partial CDS"				
							if col7 == "-":
								if dico_start[col9] in startCodons:
									start=col5
								else:
									start="<"+col5
									ext="5' Partial CDS"
								stop=col4
							col9=line.split("\t")[8].split("_")[0].rstrip()		
							list_tmp_gene.append(int(start.replace("<","").replace(">","")))
							list_tmp_gene.append(int(stop.replace("<","").replace(">","")))
							list_tmp_cds.append(start+"\t"+stop)
							dico_tmp_gff[int(col4)+1]=col1+"\t"+col2+"\t"+"CDS"+"\t"+col4+"\t"+col5+"\t"+col6+"\t"+col7+"\t"+col8+"\t"+"Parent="+col9+".1 CDS "+ext
							
					c+=1
				list_tmp_gene=sorted(list_tmp_gene)

				if col7 == "+":	
					tout.write(str(list_tmp_gene[0])+"\t"+str(list_tmp_gene[-1])+"\tgene\n")
					gout.write(col1+"\t"+col2+"\t"+"gene"+"\t"+str(list_tmp_gene[0])+"\t"+str(list_tmp_gene[-1])+"\t"+col6+"\t"+col7+"\t"+col8+"\t"+"Parent="+col9+"\n")
				if col7 == "-":	
					tout.write(str(list_tmp_gene[-1])+"\t"+str(list_tmp_gene[0])+"\tgene\n")
					gout.write(col1+"\t"+col2+"\t"+"gene"+"\t"+str(list_tmp_gene[0])+"\t"+str(list_tmp_gene[-1])+"\t"+col6+"\t"+col7+"\t"+col8+"\t"+"Parent="+col9+"\n")
				tout.write("\t\t\tgene\t"+col9+"\n")
				for line in list_tmp_cds:
					tout.write(line+"\n")

				if dico_product.has_key(col9):
					tout.write("\t\t\tproduct\t"+dico_product.get(col9)+"\n")
				else:
					tout.write("\t\t\tproduct\t"+col9+"\n")
				tout.write("\t\t\ttransl_table\t"+sys.argv[3]+"\n")
				if ext != "":
					tout.write("\t\t\tnote\t"+ext+"\n")
				sorted_y = sorted(dico_tmp_gff.items(), key=operator.itemgetter(0))
				sorted_dico = collections.OrderedDict(sorted_y)
				
				for key, val in sorted_dico.items():
					gout.write(val+"\n")
								
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
			"""if col7 == "+":
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
			gout.write(col1+"\t"+col2+"\t"+"rRNA"+"\t"+col4+"\t"+col5+"\t"+col6+"\t"+col7+"\t"+col8+"\t"+"Name="+col9+" rRNA"+"\n")"""
			if dicog.get(col9) == 1:
				if col7 == "+":
					tout.write(col4+"\t"+col5+"\t"+"gene\n")
					tout.write("\t\t\tgene\t"+col9+"\n")
					tout.write(col4+"\t"+col5+"\t"+"rRNA\n")
					tout.write("\t\t\tgene\t"+col9+"\n")
					if dico_product.has_key(col9):
						tout.write("\t\t\tproduct\t"+dico_product.get(col9)+"\n")
					else:
						tout.write("\t\t\tproduct\tunknown\n")
				if col7 == "-":
					tout.write(col5+"\t"+col4+"\t"+"gene\n")
					tout.write("\t\t\tgene\t"+col9+"\n")
					tout.write(col5+"\t"+col4+"\t"+"rRNA\n")
					tout.write("\t\t\tgene\t"+col9+"\n")
					if dico_product.has_key(col9):
						tout.write("\t\t\tproduct\t"+dico_product.get(col9)+"\n")
					else:
						tout.write("\t\t\tproduct\tunknown\n")				
				gout.write(col1+"\t"+col2+"\t"+"gene"+"\t"+col4+"\t"+col5+"\t"+col6+"\t"+col7+"\t"+col8+"\t"+"Name="+col9+" gene"+"\n")
				gout.write(col1+"\t"+col2+"\t"+"rRNA"+"\t"+col4+"\t"+col5+"\t"+col6+"\t"+col7+"\t"+col8+"\t"+"Name="+col9+" rRNA"+"\n")
			if dicog.get(col9) > 1:
				c=1
				list_tmp_gene=[]
				list_tmp_cds=[]
				dico_tmp_gff={}
				for line in dicogl.get(col9).split(";"):
					col1=seqID
					col2=line.split("\t")[1]
					col3=line.split("\t")[2]
					col4=line.split("\t")[3]
					col5=line.split("\t")[4]
					col6=line.split("\t")[5]
					col7=line.split("\t")[6]
					col8=line.split("\t")[7]
					col9=line.split("\t")[8].split("_")[0].rstrip()
					size=int(col5)-int(col4)+1
					if c == 1:
						if col7 == "+":
								start=col4
								stop=col5
						if col7 == "-":
							stop=col4
							start=col5
						list_tmp_gene.append(int(start.replace("<","").replace(">","")))
						list_tmp_gene.append(int(stop.replace("<","").replace(">","")))
						list_tmp_cds.append(start+"\t"+stop+"\t"+"rRNA")
						dico_tmp_gff[int(col4)+1]=col1+"\t"+col2+"\t"+"rRNA"+"\t"+col4+"\t"+col5+"\t"+col6+"\t"+col7+"\t"+col8+"\t"+"Parent="+col9+".1 rRNA "
							
					if c > 1 and c < int(dicog.get(col9)):
						if col7 == "+":
							start=col4
							stop=col5
						if col7 == "-":
							start=col5
							stop=col4
						col9=line.split("\t")[8].split("_")[0].rstrip()
						list_tmp_gene.append(int(start.replace("<","").replace(">","")))
						list_tmp_gene.append(int(stop.replace("<","").replace(">","")))
						list_tmp_cds.append(start+"\t"+stop)
						dico_tmp_gff[int(col4)+1]=col1+"\t"+col2+"\t"+"rRNA"+"\t"+col4+"\t"+col5+"\t"+col6+"\t"+col7+"\t"+col8+"\t"+"Parent="+col9+".1 rRNA "
							
					if c > 1 and c == int(dicog.get(col9)):
						if col7 == "+":
							start=col4
							stop=col5			
						if col7 == "-":
							start=col5
							stop=col4
						col9=line.split("\t")[8].split("_")[0].rstrip()		
						list_tmp_gene.append(int(start.replace("<","").replace(">","")))
						list_tmp_gene.append(int(stop.replace("<","").replace(">","")))
						list_tmp_cds.append(start+"\t"+stop)
						dico_tmp_gff[int(col4)+1]=col1+"\t"+col2+"\t"+"rRNA"+"\t"+col4+"\t"+col5+"\t"+col6+"\t"+col7+"\t"+col8+"\t"+"Parent="+col9+".1 rRNA "
							
					c+=1
				list_tmp_gene=sorted(list_tmp_gene)

				if col7 == "+":	
					tout.write(str(list_tmp_gene[0])+"\t"+str(list_tmp_gene[-1])+"\tgene\n")
					gout.write(col1+"\t"+col2+"\t"+"gene"+"\t"+str(list_tmp_gene[0])+"\t"+str(list_tmp_gene[-1])+"\t"+col6+"\t"+col7+"\t"+col8+"\t"+"Parent="+col9+"\n")
				if col7 == "-":	
					tout.write(str(list_tmp_gene[-1])+"\t"+str(list_tmp_gene[0])+"\tgene\n")
					gout.write(col1+"\t"+col2+"\t"+"gene"+"\t"+str(list_tmp_gene[0])+"\t"+str(list_tmp_gene[-1])+"\t"+col6+"\t"+col7+"\t"+col8+"\t"+"Parent="+col9+"\n")
				tout.write("\t\t\tgene\t"+col9+"\n")
				for line in list_tmp_cds:
					tout.write(line+"\n")

				if dico_product.has_key(col9):
					tout.write("\t\t\tproduct\t"+dico_product.get(col9)+"\n")
				else:
					tout.write("\t\t\tproduct\t"+col9+"\n")
				if ext != "":
					tout.write("\t\t\tnote\t"+ext+"\n")
				sorted_y = sorted(dico_tmp_gff.items(), key=operator.itemgetter(0))
				sorted_dico = collections.OrderedDict(sorted_y)
				
				for key, val in sorted_dico.items():
					gout.write(val+"\n")
	
		if "rrnS" in col9:
			if dicog.get(col9) == 1:
				if col7 == "+":
					tout.write(col4+"\t"+col5+"\t"+"gene\n")
					tout.write("\t\t\tgene\t"+col9+"\n")
					tout.write(col4+"\t"+col5+"\t"+"rRNA\n")
					tout.write("\t\t\tgene\t"+col9+"\n")
					if dico_product.has_key(col9):
						tout.write("\t\t\tproduct\t"+dico_product.get(col9)+"\n")
					else:
						tout.write("\t\t\tproduct\tunknown\n")
				if col7 == "-":
					tout.write(col5+"\t"+col4+"\t"+"gene\n")
					tout.write("\t\t\tgene\t"+col9+"\n")
					tout.write(col5+"\t"+col4+"\t"+"rRNA\n")
					tout.write("\t\t\tgene\t"+col9+"\n")
					if dico_product.has_key(col9):
						tout.write("\t\t\tproduct\t"+dico_product.get(col9)+"\n")
					else:
						tout.write("\t\t\tproduct\tunknown\n")				
				gout.write(col1+"\t"+col2+"\t"+"gene"+"\t"+col4+"\t"+col5+"\t"+col6+"\t"+col7+"\t"+col8+"\t"+"Name="+col9+" gene"+"\n")
				gout.write(col1+"\t"+col2+"\t"+"rRNA"+"\t"+col4+"\t"+col5+"\t"+col6+"\t"+col7+"\t"+col8+"\t"+"Name="+col9+" rRNA"+"\n")
			if dicog.get(col9) > 1:
				c=1
				list_tmp_gene=[]
				list_tmp_cds=[]
				dico_tmp_gff={}
				for line in dicogl.get(col9).split(";"):
					col1=seqID
					col2=line.split("\t")[1]
					col3=line.split("\t")[2]
					col4=line.split("\t")[3]
					col5=line.split("\t")[4]
					col6=line.split("\t")[5]
					col7=line.split("\t")[6]
					col8=line.split("\t")[7]
					col9=line.split("\t")[8].split("_")[0].rstrip()
					size=int(col5)-int(col4)+1
					if c == 1:
						if col7 == "+":
								start=col4
								stop=col5
						if col7 == "-":
							stop=col4
							start=col5
						list_tmp_gene.append(int(start.replace("<","").replace(">","")))
						list_tmp_gene.append(int(stop.replace("<","").replace(">","")))
						list_tmp_cds.append(start+"\t"+stop+"\t"+"rRNA")
						dico_tmp_gff[int(col4)+1]=col1+"\t"+col2+"\t"+"rRNA"+"\t"+col4+"\t"+col5+"\t"+col6+"\t"+col7+"\t"+col8+"\t"+"Parent="+col9+".1 rRNA "
							
					if c > 1 and c < int(dicog.get(col9)):
						if col7 == "+":
							start=col4
							stop=col5
						if col7 == "-":
							start=col5
							stop=col4
						col9=line.split("\t")[8].split("_")[0].rstrip()
						list_tmp_gene.append(int(start.replace("<","").replace(">","")))
						list_tmp_gene.append(int(stop.replace("<","").replace(">","")))
						list_tmp_cds.append(start+"\t"+stop)
						dico_tmp_gff[int(col4)+1]=col1+"\t"+col2+"\t"+"rRNA"+"\t"+col4+"\t"+col5+"\t"+col6+"\t"+col7+"\t"+col8+"\t"+"Parent="+col9+".1 rRNA "
							
					if c > 1 and c == int(dicog.get(col9)):
						if col7 == "+":
							start=col4
							stop=col5			
						if col7 == "-":
							start=col5
							stop=col4
						col9=line.split("\t")[8].split("_")[0].rstrip()		
						list_tmp_gene.append(int(start.replace("<","").replace(">","")))
						list_tmp_gene.append(int(stop.replace("<","").replace(">","")))
						list_tmp_cds.append(start+"\t"+stop)
						dico_tmp_gff[int(col4)+1]=col1+"\t"+col2+"\t"+"rRNA"+"\t"+col4+"\t"+col5+"\t"+col6+"\t"+col7+"\t"+col8+"\t"+"Parent="+col9+".1 rRNA "
							
					c+=1
				list_tmp_gene=sorted(list_tmp_gene)

				if col7 == "+":	
					tout.write(str(list_tmp_gene[0])+"\t"+str(list_tmp_gene[-1])+"\tgene\n")
					gout.write(col1+"\t"+col2+"\t"+"gene"+"\t"+str(list_tmp_gene[0])+"\t"+str(list_tmp_gene[-1])+"\t"+col6+"\t"+col7+"\t"+col8+"\t"+"Parent="+col9+"\n")
				if col7 == "-":	
					tout.write(str(list_tmp_gene[-1])+"\t"+str(list_tmp_gene[0])+"\tgene\n")
					gout.write(col1+"\t"+col2+"\t"+"gene"+"\t"+str(list_tmp_gene[0])+"\t"+str(list_tmp_gene[-1])+"\t"+col6+"\t"+col7+"\t"+col8+"\t"+"Parent="+col9+"\n")
				tout.write("\t\t\tgene\t"+col9+"\n")
				for line in list_tmp_cds:
					tout.write(line+"\n")

				if dico_product.has_key(col9):
					tout.write("\t\t\tproduct\t"+dico_product.get(col9)+"\n")
				else:
					tout.write("\t\t\tproduct\t"+col9+"\n")
				if ext != "":
					tout.write("\t\t\tnote\t"+ext+"\n")
				sorted_y = sorted(dico_tmp_gff.items(), key=operator.itemgetter(0))
				sorted_dico = collections.OrderedDict(sorted_y)
				
				for key, val in sorted_dico.items():
					gout.write(val+"\n")
gout.close()
tout.close()

	
	
	
