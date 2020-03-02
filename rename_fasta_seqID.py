#!/usr/bin/python
#Version: 1.2
#Authors: Allio Remi & Schomaker-Bastos Alex
#ISEM - CNRS - LAMPADA - IBQM - UFRJ

import sys
import os.path

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
fasta = open(sys.argv[2])
fout = open(sys.argv[3], 'w')
c=sys.argv[4]

for name, seq in read_fasta(fasta):
	fout.write(">"+seqID+"."+str(c)+"\n"+seq+"\n")
