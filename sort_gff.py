#! /usr/bin/env python

import sys
import os.path
import operator
import collections

file1=open(sys.argv[1])

dico={}
dicotrna={}
dicof={}
c=0

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
fout=sys.argv[1].split("_raw.gff")[0]+"_filtered.gff"
fout=open(fout, "w")

for k, v in sorted_dict.items():
	if not dicotrna.has_key(v.split("\t")[8]):
		fout.write(v+"\n")
	
	
	
	
	
