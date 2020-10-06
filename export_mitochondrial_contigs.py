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
import sys
import os.path


f=open(sys.argv[1])
g=open(sys.argv[1])
evalmin=sys.argv[2]

# check species are present only once
dico={}
dicof={}
for line in f:
	trans=line.split("	")[1]
	trans=trans.split("	")[0]
	valeur=trans
	dico[trans]=0
	dicof[valeur]=str("0 ()")



for cle, valeur in dico.items():
	for line in open(sys.argv[1]):
		if cle in line:
			blast=line.split("	")[2]
			blast=blast.split("	")[0]
			cibl=line.split("	")[0]
			blast=float(blast)
			eVal=line.split("	")[10]
			eVal=eVal.split("	")[0]
			test=dicof.get(cle)
			test=test.split(" ")[0]
			if blast > float(test) and float(eVal) > float(evalmin):
				dicof[cle]=str(float(blast))+" ("+cibl+")"
for cle, valeur in dicof.items():
	print cle+" === "+str(valeur)
	


