#!/usr/bin/python
#Version: 1.2
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
import argparse, os, shlex, shutil, sys
import numpy as np
import pandas as pd
import glob

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

class SmartFormatter(argparse.HelpFormatter):

    def _split_lines(self, text, width):
        if text.startswith('R|'):
            return text[2:].splitlines()  
        # this is the RawTextHelpFormatter._split_lines
        return argparse.HelpFormatter._split_lines(self, text, width)

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='This script is a subpart of Mitofinder designed to create files for tbl2asn submission.', formatter_class=SmartFormatter)
	parser.add_argument('-i', '--index', help='Index file (comma-delimited table)', default="", dest='table')
	parser.add_argument('-v', '--version', help="MitoFinder version 1.2", default=False, dest='versionCheck', action='store_true')
	
	args = parser.parse_args()

	if args.versionCheck == True:
		print "MitoFinder version 1.2"
		exit()


	table = pd.read_csv(args.table)
	df=pd.DataFrame(table)
	size=len(list(df.columns))
	dico={}
	initial_path=os.getcwd()+"/"
	for line in list(df.columns):
		dico[line]=line
	
	if not dico.has_key("Directory path"):
		print "ERROR: The column \"Directory path\" is not found in the index files"
		exit()
	if not dico.has_key("Seq ID"):
		print "ERROR: The column \"Seq ID\" is not found in the index files"
		exit()
	
	dico_seqid={}
	dico_sra1={}
	dico_sra2={}
	
	for index, row in df.iterrows():
		seqid=str(row["Seq ID"])
		dico_seqid[seqid]=""
		dico_sra1[seqid]=""
		dico_sra2[seqid]=""
		for col, value in sorted(dico.items()):
			if col != "Directory path" and col != "Seq ID":
				if not row[col] != row[col]:
					dico_seqid[seqid]=str(dico_seqid.get(seqid))+" ["+str(col)+"="+str(row[col])+"]"		
							
	for index, row in df.iterrows():
		path=row["Directory path"]
		seqid=str(row["Seq ID"])
		if path != path:
			path=os.getcwd()+"/"
		elif path == "NA":
			path=os.getcwd()+"/"
		else:
			path=os.path.join(initial_path,os.path.expanduser(path))
		seq=os.path.join(path,seqid)
		fasta=os.path.join(path,seqid+".fasta")
		if os.path.isfile(os.path.join(path,seqid+".fsa")):
			os.remove(os.path.join(path,seqid+".fsa"))
		if os.path.isfile(os.path.join(path,seqid+".tbl")):
			os.remove(os.path.join(path,seqid+".tbl"))
		c=0
		print "Formatting "+str(seqid)+" FASTA file\n"
		for f in sorted(glob.glob(seq+"*.fasta")):
			c += 1
			tmp=open(f)
			fout = os.path.join(path,seqid+".fsa")
			fout=open(fout,"a")
			for name, seq in read_fasta(tmp):
				name=str(name)+str(dico_seqid.get(seqid))
				fout.write(name+"\n"+seq+"\n")	
			if os.path.isfile(f.split(".fasta")[0]+".tbl"):
				tblout=os.path.join(path,seqid+".tbl")
				tblout=open(tblout,"a")
				tmp=open(f.split(".fasta")[0]+".tbl")
				for line in tmp:
					tblout.write(line.rstrip()+"\n")
				tmp.close()
				tblout.close()
			else:
				print "WARNING: no tbl files (.tbl) were found for "+str(seqid)+ " at the following path:\n"+str(path)
				print "The annoation inforamtions for the contig \""+str(f.split("/")[-1])+ "\" will not be found in the final tbl file : \""+ str(f.split(".fasta")[0].split("/")[-1]+".tbl")+"\"\n"
			fout.close()
		if c == 0:
			print "WARNING: no fasta files (.fasta) were found for "+str(seqid)+ " at the following path:\n"+str(path)+"\n"
			print "WARNING: no tbl files (.tbl) were found for "+str(seqid)+ " at the following path:\n"+str(path)+"\n"













