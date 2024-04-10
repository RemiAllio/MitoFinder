#!/usr/bin/python2.7
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

import argparse, os, shlex, shutil, sys
import geneChecker, genbankOutput, runMegahit,circularizationCheck
import subprocess
from subprocess import Popen
from Bio import SeqIO, SeqFeature, SeqUtils
from Bio.Alphabet import generic_dna, generic_protein
import glob


record=SeqIO.read(sys.argv[1], "genbank", generic_dna)
ID=sys.argv[3]

with open(sys.argv[2], 'a') as importantFeaturesFile:
    for feature in record.features:
        if feature.type.lower() == 'cds':
            if 'gene' in feature.qualifiers:
                featureName = feature.qualifiers['gene'][0]
            elif 'product' in feature.qualifiers:
                featureName = feature.qualifiers['product'][0]
                featureName = ''.join(featureName.split())
                                        
            importantFeaturesFile.write('>' + str(ID) + "@" + featureName + '\n')
            importantFeaturesFile.write(str(feature.extract(record).seq)+'\n')

        if feature.type.lower() == 'rrna':
            if 'gene' in feature.qualifiers:
                featureName = feature.qualifiers['gene'][0]
            elif 'product' in feature.qualifiers:
                featureName = feature.qualifiers['product'][0]
                featureName = ''.join(featureName.split())
            
            importantFeaturesFile.write('>' + str(ID) + "@" + featureName + '\n')
            importantFeaturesFile.write(str(feature.extract(record).seq)+'\n')
