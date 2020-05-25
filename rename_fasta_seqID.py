#!/usr/bin/python
#Version: 1.2
#Authors: Allio Remi & Schomaker-Bastos Alex
#ISEM - CNRS - LAMPADA - IBQM - UFRJ

import sys
import os.path
from Bio import SeqIO, SeqFeature
from Bio.Alphabet import generic_dna, generic_protein, IUPAC
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

seqID=sys.argv[1]
fout = open(sys.argv[3], 'w')
c=sys.argv[4]
direction=sys.argv[5]

resultFile=SeqIO.read(open(sys.argv[2], 'rU'), "fasta", generic_dna)
if direction == "+":
	if sys.argv[6] == "yes":
		fout.write(">"+seqID+"."+str(c)+"\n"+str(resultFile.seq)+"\n")
	else:
		fout.write(">"+str(resultFile.id)+"\n"+str(resultFile.seq)+"\n")
else:
	if sys.argv[6] == "yes":
		fout.write(">"+seqID+"."+str(c)+"\n"+str(resultFile.seq.reverse_complement())+"\n")
	else:
		fout.write(">"+str(resultFile.id)+" (reverse)\n"+str(resultFile.seq)+"\n")
