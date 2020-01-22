#! /usr/bin/env python

import sys
import os.path
import operator
import collections

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
gout=sys.argv[1].split("_raw.gff")[0]+"_final.gff"
gout=open(gout, "w")
gout.write(seqID+"\t"+"mitofinder\tsource\t1\t"+str(length)+"\t.\t+\t.\tName=source "+seqID.replace("_"," ")+"\n")
tout=sys.argv[1].split("_raw.gff")[0]+"_final.tbl"
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
		if "COX1" == col9:
			size=int(col5)-int(col4)+1
			if size%3 == 0:
				ext=""
				tout.write(col4+"\t"+col5+"\t"+"gene\n")
				tout.write("\t\t\tgene\t"+col9+"\n")
				tout.write(col4+"\t"+col5+"\t"+"CDS\n")
				tout.write("\t\t\tproduct\tcytochrome c oxidase subunit I\n")
				tout.write("\t\t\ttransl_table\t"+sys.argv[3]+"\n")
			else:
				ext=" Note: Not a multiple of 3"
				if col7 == "+":
					tout.write(col4+"\t"+col5+">"+"\t"+"gene\n")
					tout.write("\t\t\tgene\t"+col9+"\n")
					tout.write(col4+"\t"+col5+">"+"\t"+"CDS\n")
					tout.write("\t\t\tproduct\tcytochrome c oxidase subunit I\n")
					tout.write("\t\t\ttransl_table\t"+sys.argv[3]+"\n")
					tout.write("\t\t\tNote\tPartial sequence\n")
				if col7 == "-":
					tout.write("<"+col4+"\t"+col5+"\t"+"gene\n")
					tout.write("\t\t\tgene\t"+col9+"\n")
					tout.write("<"+col4+"\t"+col5+"\t"+"CDS\n")
					tout.write("\t\t\tproduct\tcytochrome c oxidase subunit I\n")
					tout.write("\t\t\ttransl_table\t"+sys.argv[3]+"\n")
					tout.write("\t\t\tNote\tPartial sequence\n")
			gout.write(col1+"\t"+col2+"\t"+"gene"+"\t"+col4+"\t"+col5+"\t"+col6+"\t"+col7+"\t"+col8+"\t"+"Name="+col9+" gene"+ext+"\n")
			gout.write(col1+"\t"+col2+"\t"+"CDS"+"\t"+col4+"\t"+col5+"\t"+col6+"\t"+col7+"\t"+col8+"\t"+"Name="+col9+" CDS"+ext+"\n")
		
		if "COX2" == col9:
			size=int(col5)-int(col4)+1
			if size%3 == 0:
				ext=""
				tout.write(col4+"\t"+col5+"\t"+"gene\n")
				tout.write("\t\t\tgene\t"+col9+"\n")
				tout.write(col4+"\t"+col5+"\t"+"CDS\n")
				tout.write("\t\t\tproduct\tcytochrome c oxidase subunit II\n")
				tout.write("\t\t\ttransl_table\t"+sys.argv[3]+"\n")
			else:
				ext=" Note: Not a multiple of 3"
				if col7 == "+":
					tout.write(col4+"\t"+col5+">"+"\t"+"gene\n")
					tout.write("\t\t\tgene\t"+col9+"\n")
					tout.write(col4+"\t"+col5+">"+"\t"+"CDS\n")
					tout.write("\t\t\tproduct\tcytochrome c oxidase subunit II\n")
					tout.write("\t\t\ttransl_table\t"+sys.argv[3]+"\n")
					tout.write("\t\t\tNote\tPartial sequence\n")
				if col7 == "-":
					tout.write("<"+col4+"\t"+col5+"\t"+"gene\n")
					tout.write("\t\t\tgene\t"+col9+"\n")
					tout.write("<"+col4+"\t"+col5+"\t"+"CDS\n")
					tout.write("\t\t\tproduct\tcytochrome c oxidase subunit II\n")
					tout.write("\t\t\ttransl_table\t"+sys.argv[3]+"\n")
					tout.write("\t\t\tNote\tPartial sequence\n")
			gout.write(col1+"\t"+col2+"\t"+"gene"+"\t"+col4+"\t"+col5+"\t"+col6+"\t"+col7+"\t"+col8+"\t"+"Name="+col9+" gene"+ext+"\n")
			gout.write(col1+"\t"+col2+"\t"+"CDS"+"\t"+col4+"\t"+col5+"\t"+col6+"\t"+col7+"\t"+col8+"\t"+"Name="+col9+" CDS"+ext+"\n")
		
		if "COX3" == col9:
			size=int(col5)-int(col4)+1
			if size%3 == 0:
				ext=""
				tout.write(col4+"\t"+col5+"\t"+"gene\n")
				tout.write("\t\t\tgene\t"+col9+"\n")
				tout.write(col4+"\t"+col5+"\t"+"CDS\n")
				tout.write("\t\t\tproduct\tcytochrome c oxidase subunit III\n")
				tout.write("\t\t\ttransl_table\t"+sys.argv[3]+"\n")
			else:
				ext=" Note: Not a multiple of 3"
				if col7 == "+":
					tout.write(col4+"\t"+col5+">"+"\t"+"gene\n")
					tout.write("\t\t\tgene\t"+col9+"\n")
					tout.write(col4+"\t"+col5+">"+"\t"+"CDS\n")
					tout.write("\t\t\tproduct\tcytochrome c oxidase subunit III\n")
					tout.write("\t\t\ttransl_table\t"+sys.argv[3]+"\n")
					tout.write("\t\t\tNote\tPartial sequence\n")
				if col7 == "-":
					tout.write("<"+col4+"\t"+col5+"\t"+"gene\n")
					tout.write("\t\t\tgene\t"+col9+"\n")
					tout.write("<"+col4+"\t"+col5+"\t"+"CDS\n")
					tout.write("\t\t\tproduct\tcytochrome c oxidase subunit III\n")
					tout.write("\t\t\ttransl_table\t"+sys.argv[3]+"\n")
					tout.write("\t\t\tNote\tPartial sequence\n")
			gout.write(col1+"\t"+col2+"\t"+"gene"+"\t"+col4+"\t"+col5+"\t"+col6+"\t"+col7+"\t"+col8+"\t"+"Name="+col9+" gene"+ext+"\n")
			gout.write(col1+"\t"+col2+"\t"+"CDS"+"\t"+col4+"\t"+col5+"\t"+col6+"\t"+col7+"\t"+col8+"\t"+"Name="+col9+" CDS"+ext+"\n")
	
		if "ND1" == col9:
			size=int(col5)-int(col4)+1
			if size%3 == 0:
				ext=""
				tout.write(col4+"\t"+col5+"\t"+"gene\n")
				tout.write("\t\t\tgene\t"+col9+"\n")
				tout.write(col4+"\t"+col5+"\t"+"CDS\n")
				tout.write("\t\t\tproduct\tNADH dehydrogenase subunit 1\n")
				tout.write("\t\t\ttransl_table\t"+sys.argv[3]+"\n")
			else:
				ext=" Note: Not a multiple of 3"
				if col7 == "+":
					tout.write(col4+"\t"+col5+">"+"\t"+"gene\n")
					tout.write("\t\t\tgene\t"+col9+"\n")
					tout.write(col4+"\t"+col5+">"+"\t"+"CDS\n")
					tout.write("\t\t\tproduct\tNADH dehydrogenase subunit 1\n")
					tout.write("\t\t\ttransl_table\t"+sys.argv[3]+"\n")
					tout.write("\t\t\tNote\tPartial sequence\n")
				if col7 == "-":
					tout.write("<"+col4+"\t"+col5+"\t"+"gene\n")
					tout.write("\t\t\tgene\t"+col9+"\n")
					tout.write("<"+col4+"\t"+col5+"\t"+"CDS\n")
					tout.write("\t\t\tproduct\tNADH dehydrogenase subunit 1\n")
					tout.write("\t\t\ttransl_table\t"+sys.argv[3]+"\n")
					tout.write("\t\t\tNote\tPartial sequence\n")
			gout.write(col1+"\t"+col2+"\t"+"gene"+"\t"+col4+"\t"+col5+"\t"+col6+"\t"+col7+"\t"+col8+"\t"+"Name="+col9+" gene"+ext+"\n")
			gout.write(col1+"\t"+col2+"\t"+"CDS"+"\t"+col4+"\t"+col5+"\t"+col6+"\t"+col7+"\t"+col8+"\t"+"Name="+col9+" CDS"+ext+"\n")

		if "ND2" == col9:
			size=int(col5)-int(col4)+1
			if size%3 == 0:
				ext=""
				tout.write(col4+"\t"+col5+"\t"+"gene\n")
				tout.write("\t\t\tgene\t"+col9+"\n")
				tout.write(col4+"\t"+col5+"\t"+"CDS\n")
				tout.write("\t\t\tproduct\tNADH dehydrogenase subunit 2\n")
				tout.write("\t\t\ttransl_table\t"+sys.argv[3]+"\n")
			else:
				ext=" Note: Not a multiple of 3"
				if col7 == "+":
					tout.write(col4+"\t"+col5+">"+"\t"+"gene\n")
					tout.write("\t\t\tgene\t"+col9+"\n")
					tout.write(col4+"\t"+col5+">"+"\t"+"CDS\n")
					tout.write("\t\t\tproduct\tNADH dehydrogenase subunit 2\n")
					tout.write("\t\t\ttransl_table\t"+sys.argv[3]+"\n")
					tout.write("\t\t\tNote\tPartial sequence\n")
				if col7 == "-":
					tout.write("<"+col4+"\t"+col5+"\t"+"gene\n")
					tout.write("\t\t\tgene\t"+col9+"\n")
					tout.write("<"+col4+"\t"+col5+"\t"+"CDS\n")
					tout.write("\t\t\tproduct\tNADH dehydrogenase subunit 2\n")
					tout.write("\t\t\ttransl_table\t"+sys.argv[3]+"\n")
					tout.write("\t\t\tNote\tPartial sequence\n")
			gout.write(col1+"\t"+col2+"\t"+"gene"+"\t"+col4+"\t"+col5+"\t"+col6+"\t"+col7+"\t"+col8+"\t"+"Name="+col9+" gene"+ext+"\n")
			gout.write(col1+"\t"+col2+"\t"+"CDS"+"\t"+col4+"\t"+col5+"\t"+col6+"\t"+col7+"\t"+col8+"\t"+"Name="+col9+" CDS"+ext+"\n")
		
		if "ND3" == col9:
			size=int(col5)-int(col4)+1
			if size%3 == 0:
				ext=""
				tout.write(col4+"\t"+col5+"\t"+"gene\n")
				tout.write("\t\t\tgene\t"+col9+"\n")
				tout.write(col4+"\t"+col5+"\t"+"CDS\n")
				tout.write("\t\t\tproduct\tNADH dehydrogenase subunit 3\n")
				tout.write("\t\t\ttransl_table\t"+sys.argv[3]+"\n")
			else:
				ext=" Note: Not a multiple of 3"
				if col7 == "+":
					tout.write(col4+"\t"+col5+">"+"\t"+"gene\n")
					tout.write("\t\t\tgene\t"+col9+"\n")
					tout.write(col4+"\t"+col5+">"+"\t"+"CDS\n")
					tout.write("\t\t\tproduct\tNADH dehydrogenase subunit 3\n")
					tout.write("\t\t\ttransl_table\t"+sys.argv[3]+"\n")
					tout.write("\t\t\tNote\tPartial sequence\n")
				if col7 == "-":
					tout.write("<"+col4+"\t"+col5+"\t"+"gene\n")
					tout.write("\t\t\tgene\t"+col9+"\n")
					tout.write("<"+col4+"\t"+col5+"\t"+"CDS\n")
					tout.write("\t\t\tproduct\tNADH dehydrogenase subunit 3\n")
					tout.write("\t\t\ttransl_table\t"+sys.argv[3]+"\n")
					tout.write("\t\t\tNote\tPartial sequence\n")
			gout.write(col1+"\t"+col2+"\t"+"gene"+"\t"+col4+"\t"+col5+"\t"+col6+"\t"+col7+"\t"+col8+"\t"+"Name="+col9+" gene"+ext+"\n")
			gout.write(col1+"\t"+col2+"\t"+"CDS"+"\t"+col4+"\t"+col5+"\t"+col6+"\t"+col7+"\t"+col8+"\t"+"Name="+col9+" CDS"+ext+"\n")

		if "ND4" == col9:
			size=int(col5)-int(col4)+1
			if size%3 == 0:
				ext=""
				tout.write(col4+"\t"+col5+"\t"+"gene\n")
				tout.write("\t\t\tgene\t"+col9+"\n")
				tout.write(col4+"\t"+col5+"\t"+"CDS\n")
				tout.write("\t\t\tproduct\tNADH dehydrogenase subunit 4\n")
				tout.write("\t\t\ttransl_table\t"+sys.argv[3]+"\n")
			else:
				ext=" Note: Not a multiple of 3"
				if col7 == "+":
					tout.write(col4+"\t"+col5+">"+"\t"+"gene\n")
					tout.write("\t\t\tgene\t"+col9+"\n")
					tout.write(col4+"\t"+col5+">"+"\t"+"CDS\n")
					tout.write("\t\t\tproduct\tNADH dehydrogenase subunit 4\n")
					tout.write("\t\t\ttransl_table\t"+sys.argv[3]+"\n")
					tout.write("\t\t\tNote\tPartial sequence\n")
				if col7 == "-":
					tout.write("<"+col4+"\t"+col5+"\t"+"gene\n")
					tout.write("\t\t\tgene\t"+col9+"\n")
					tout.write("<"+col4+"\t"+col5+"\t"+"CDS\n")
					tout.write("\t\t\tproduct\tNADH dehydrogenase subunit 4\n")
					tout.write("\t\t\ttransl_table\t"+sys.argv[3]+"\n")
					tout.write("\t\t\tNote\tPartial sequence\n")
			gout.write(col1+"\t"+col2+"\t"+"gene"+"\t"+col4+"\t"+col5+"\t"+col6+"\t"+col7+"\t"+col8+"\t"+"Name="+col9+" gene"+ext+"\n")
			gout.write(col1+"\t"+col2+"\t"+"CDS"+"\t"+col4+"\t"+col5+"\t"+col6+"\t"+col7+"\t"+col8+"\t"+"Name="+col9+" CDS"+ext+"\n")

		if "ND4L" == col9:
			size=int(col5)-int(col4)+1
			if size%3 == 0:
				ext=""
				tout.write(col4+"\t"+col5+"\t"+"gene\n")
				tout.write("\t\t\tgene\t"+col9+"\n")
				tout.write(col4+"\t"+col5+"\t"+"CDS\n")
				tout.write("\t\t\tproduct\tNADH dehydrogenase subunit 4L\n")
				tout.write("\t\t\ttransl_table\t"+sys.argv[3]+"\n")
			else:
				ext=" Note: Not a multiple of 3"
				if col7 == "+":
					tout.write(col4+"\t"+col5+">"+"\t"+"gene\n")
					tout.write("\t\t\tgene\t"+col9+"\n")
					tout.write(col4+"\t"+col5+">"+"\t"+"CDS\n")
					tout.write("\t\t\tproduct\tNADH dehydrogenase subunit 4L\n")
					tout.write("\t\t\ttransl_table\t"+sys.argv[3]+"\n")
					tout.write("\t\t\tNote\tPartial sequence\n")
				if col7 == "-":
					tout.write("<"+col4+"\t"+col5+"\t"+"gene\n")
					tout.write("\t\t\tgene\t"+col9+"\n")
					tout.write("<"+col4+"\t"+col5+"\t"+"CDS\n")
					tout.write("\t\t\tproduct\tNADH dehydrogenase subunit 4L\n")
					tout.write("\t\t\ttransl_table\t"+sys.argv[3]+"\n")
					tout.write("\t\t\tNote\tPartial sequence\n")
			gout.write(col1+"\t"+col2+"\t"+"gene"+"\t"+col4+"\t"+col5+"\t"+col6+"\t"+col7+"\t"+col8+"\t"+"Name="+col9+" gene"+ext+"\n")
			gout.write(col1+"\t"+col2+"\t"+"CDS"+"\t"+col4+"\t"+col5+"\t"+col6+"\t"+col7+"\t"+col8+"\t"+"Name="+col9+" CDS"+ext+"\n")

		if "ND5" == col9:
			size=int(col5)-int(col4)+1
			if size%3 == 0:
				ext=""
				tout.write(col4+"\t"+col5+"\t"+"gene\n")
				tout.write("\t\t\tgene\t"+col9+"\n")
				tout.write(col4+"\t"+col5+"\t"+"CDS\n")
				tout.write("\t\t\tproduct\tNADH dehydrogenase subunit 5\n")
				tout.write("\t\t\ttransl_table\t"+sys.argv[3]+"\n")
			else:
				ext=" Note: Not a multiple of 3"
				if col7 == "+":
					tout.write(col4+"\t"+col5+">"+"\t"+"gene\n")
					tout.write("\t\t\tgene\t"+col9+"\n")
					tout.write(col4+"\t"+col5+">"+"\t"+"CDS\n")
					tout.write("\t\t\tproduct\tNADH dehydrogenase subunit 5\n")
					tout.write("\t\t\ttransl_table\t"+sys.argv[3]+"\n")
					tout.write("\t\t\tNote\tPartial sequence\n")
				if col7 == "-":
					tout.write("<"+col4+"\t"+col5+"\t"+"gene\n")
					tout.write("\t\t\tgene\t"+col9+"\n")
					tout.write("<"+col4+"\t"+col5+"\t"+"CDS\n")
					tout.write("\t\t\tproduct\tNADH dehydrogenase subunit 5\n")
					tout.write("\t\t\ttransl_table\t"+sys.argv[3]+"\n")
					tout.write("\t\t\tNote\tPartial sequence\n")
			gout.write(col1+"\t"+col2+"\t"+"gene"+"\t"+col4+"\t"+col5+"\t"+col6+"\t"+col7+"\t"+col8+"\t"+"Name="+col9+" gene"+ext+"\n")
			gout.write(col1+"\t"+col2+"\t"+"CDS"+"\t"+col4+"\t"+col5+"\t"+col6+"\t"+col7+"\t"+col8+"\t"+"Name="+col9+" CDS"+ext+"\n")

		if "CYTB" == col9:
			size=int(col5)-int(col4)+1
			if size%3 == 0:
				ext=""
				tout.write(col4+"\t"+col5+"\t"+"gene\n")
				tout.write("\t\t\tgene\t"+col9+"\n")
				tout.write(col4+"\t"+col5+"\t"+"CDS\n")
				tout.write("\t\t\tproduct\tcytochrome b\n")
				tout.write("\t\t\ttransl_table\t"+sys.argv[3]+"\n")
			else:
				ext=" Note: Not a multiple of 3"
				if col7 == "+":
					tout.write(col4+"\t"+col5+">"+"\t"+"gene\n")
					tout.write("\t\t\tgene\t"+col9+"\n")
					tout.write(col4+"\t"+col5+">"+"\t"+"CDS\n")
					tout.write("\t\t\tproduct\tcytochrome b\n")
					tout.write("\t\t\ttransl_table\t"+sys.argv[3]+"\n")
					tout.write("\t\t\tNote\tPartial sequence\n")
				if col7 == "-":
					tout.write("<"+col4+"\t"+col5+"\t"+"gene\n")
					tout.write("\t\t\tgene\t"+col9+"\n")
					tout.write("<"+col4+"\t"+col5+"\t"+"CDS\n")
					tout.write("\t\t\tproduct\tcytochrome b\n")
					tout.write("\t\t\ttransl_table\t"+sys.argv[3]+"\n")
					tout.write("\t\t\tNote\tPartial sequence\n")
			gout.write(col1+"\t"+col2+"\t"+"gene"+"\t"+col4+"\t"+col5+"\t"+col6+"\t"+col7+"\t"+col8+"\t"+"Name="+col9+" gene"+ext+"\n")
			gout.write(col1+"\t"+col2+"\t"+"CDS"+"\t"+col4+"\t"+col5+"\t"+col6+"\t"+col7+"\t"+col8+"\t"+"Name="+col9+" CDS"+ext+"\n")

		if "ATP8" == col9:
			size=int(col5)-int(col4)+1
			if size%3 == 0:
				ext=""
				tout.write(col4+"\t"+col5+"\t"+"gene\n")
				tout.write("\t\t\tgene\t"+col9+"\n")
				tout.write(col4+"\t"+col5+"\t"+"CDS\n")
				tout.write("\t\t\tproduct\tATP synthase F0 subunit 8\n")
				tout.write("\t\t\ttransl_table\t"+sys.argv[3]+"\n")
			else:
				ext=" Note: Not a multiple of 3"
				if col7 == "+":
					tout.write(col4+"\t"+col5+">"+"\t"+"gene\n")
					tout.write("\t\t\tgene\t"+col9+"\n")
					tout.write(col4+"\t"+col5+">"+"\t"+"CDS\n")
					tout.write("\t\t\tproduct\tATP synthase F0 subunit 8\n")
					tout.write("\t\t\ttransl_table\t"+sys.argv[3]+"\n")
					tout.write("\t\t\tNote\tPartial sequence\n")
				if col7 == "-":
					tout.write("<"+col4+"\t"+col5+"\t"+"gene\n")
					tout.write("\t\t\tgene\t"+col9+"\n")
					tout.write("<"+col4+"\t"+col5+"\t"+"CDS\n")
					tout.write("\t\t\tproduct\tATP synthase F0 subunit 8\n")
					tout.write("\t\t\ttransl_table\t"+sys.argv[3]+"\n")
					tout.write("\t\t\tNote\tPartial sequence\n")
			gout.write(col1+"\t"+col2+"\t"+"gene"+"\t"+col4+"\t"+col5+"\t"+col6+"\t"+col7+"\t"+col8+"\t"+"Name="+col9+" gene"+ext+"\n")
			gout.write(col1+"\t"+col2+"\t"+"CDS"+"\t"+col4+"\t"+col5+"\t"+col6+"\t"+col7+"\t"+col8+"\t"+"Name="+col9+" CDS"+ext+"\n")

		if "ATP6" == col9:
			size=int(col5)-int(col4)+1
			if size%3 == 0:
				ext=""
				tout.write(col4+"\t"+col5+"\t"+"gene\n")
				tout.write("\t\t\tgene\t"+col9+"\n")
				tout.write(col4+"\t"+col5+"\t"+"CDS\n")
				tout.write("\t\t\tproduct\tATP synthase F0 subunit 6\n")
				tout.write("\t\t\ttransl_table\t"+sys.argv[3]+"\n")
			else:
				ext=" Note: Not a multiple of 3"
				if col7 == "+":
					tout.write(col4+"\t"+col5+">"+"\t"+"gene\n")
					tout.write("\t\t\tgene\t"+col9+"\n")
					tout.write(col4+"\t"+col5+">"+"\t"+"CDS\n")
					tout.write("\t\t\tproduct\tATP synthase F0 subunit 6\n")
					tout.write("\t\t\ttransl_table\t"+sys.argv[3]+"\n")
					tout.write("\t\t\tNote\tPartial sequence\n")
				if col7 == "-":
					tout.write("<"+col4+"\t"+col5+"\t"+"gene\n")
					tout.write("\t\t\tgene\t"+col9+"\n")
					tout.write("<"+col4+"\t"+col5+"\t"+"CDS\n")
					tout.write("\t\t\tproduct\tATP synthase F0 subunit 6\n")
					tout.write("\t\t\ttransl_table\t"+sys.argv[3]+"\n")
					tout.write("\t\t\tNote\tPartial sequence\n")
			gout.write(col1+"\t"+col2+"\t"+"gene"+"\t"+col4+"\t"+col5+"\t"+col6+"\t"+col7+"\t"+col8+"\t"+"Name="+col9+" gene"+ext+"\n")
			gout.write(col1+"\t"+col2+"\t"+"CDS"+"\t"+col4+"\t"+col5+"\t"+col6+"\t"+col7+"\t"+col8+"\t"+"Name="+col9+" CDS"+ext+"\n")

		if "tRNA" in col9:
			size=int(col5)-int(col4)+1
			if size <= 120:
				tout.write(col4+"\t"+col5+"\t"+"gene\n")
				tout.write("\t\t\tgene\t"+col9+"\n")
				tout.write(col4+"\t"+col5+"\t"+"tRNA\n")
				tout.write("\t\t\ttRNA\t"+col9+"\n")
				gout.write(col1+"\t"+col2+"\t"+"gene"+"\t"+col4+"\t"+col5+"\t"+col6+"\t"+col7+"\t"+col8+"\t"+"Name="+col9+" gene"+"\n")
				gout.write(col1+"\t"+col2+"\t"+"tRNA"+"\t"+col4+"\t"+col5+"\t"+col6+"\t"+col7+"\t"+col8+"\t"+"Name="+col9+" "+"\n")
		if "trn" in col9:
			size=int(col5)-int(col4)+1
			if size <= 120:
				tout.write(col4+"\t"+col5+"\t"+"gene\n")
				tout.write("\t\t\tgene\t"+col9+"\n")
				tout.write(col4+"\t"+col5+"\t"+"tRNA\n")
				tout.write("\t\t\ttRNA\t"+col9+"\n")
				gout.write(col1+"\t"+col2+"\t"+"gene"+"\t"+col4+"\t"+col5+"\t"+col6+"\t"+col7+"\t"+col8+"\t"+"Name="+col9+" gene"+"\n")
				gout.write(col1+"\t"+col2+"\t"+"tRNA"+"\t"+col4+"\t"+col5+"\t"+col6+"\t"+col7+"\t"+col8+"\t"+"Name="+col9+" tRNA"+"\n")
		if "rrnL" in col9:
			tout.write(col4+"\t"+col5+"\t"+"gene\n")
			tout.write("\t\t\tgene\t"+col9+"\n")
			tout.write(col4+"\t"+col5+"\t"+"rRNA\n")
			tout.write("\t\t\trRNA\t"+col9+"\n")
			tout.write("\t\t\tproduct\tlarge subunit ribosomal RNA\n")
			gout.write(col1+"\t"+col2+"\t"+"gene"+"\t"+col4+"\t"+col5+"\t"+col6+"\t"+col7+"\t"+col8+"\t"+"Name="+col9+" gene"+"\n")
			gout.write(col1+"\t"+col2+"\t"+"rRNA"+"\t"+col4+"\t"+col5+"\t"+col6+"\t"+col7+"\t"+col8+"\t"+"Name="+col9+" rRNA"+"\n")
		if "rrnS" in col9:
			tout.write(col4+"\t"+col5+"\t"+"gene\n")
			tout.write("\t\t\tgene\t"+col9+"\n")
			tout.write(col4+"\t"+col5+"\t"+"rRNA\n")
			tout.write("\t\t\trRNA\t"+col9+"\n")
			tout.write("\t\t\tproduct\tsmall subunit ribosomal RNA\n")
			gout.write(col1+"\t"+col2+"\t"+"gene"+"\t"+col4+"\t"+col5+"\t"+col6+"\t"+col7+"\t"+col8+"\t"+"Name="+col9+" gene"+"\n")
			gout.write(col1+"\t"+col2+"\t"+"rRNA"+"\t"+col4+"\t"+col5+"\t"+col6+"\t"+col7+"\t"+col8+"\t"+"Name="+col9+" rRNA"+"\n")
gout.close()
tout.close()

	
	
	
