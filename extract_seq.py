#!/usr/bin/python
#Version: 1.3
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

import sys
import os.path
import shlex, subprocess

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

allCDS=open(sys.argv[1])

dico={}
for name, seq in read_fasta(allCDS):
    gene=name.split("@")[1]
    sp=name.split("@")[0]
    if dico.has_key(gene):
        if len(seq) > len(dico.get(gene)):
            dico[gene]=seq
    else:
        dico[gene]=seq

for cle, valeur in dico.items():
    fout=open(cle+"_all_sp.fasta","a")
    fout.write(sp+"\n"+valeur+"\n")
