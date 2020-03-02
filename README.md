# MitoFinder version 1.2		
Allio, R., Schomaker-Bastos, A., Romiguier, J., Prosdocimi, F., Nabholz, B., & Delsuc, F.

<p align="center">
  <img src="/image/logo.png" alt="Drawing" width="250"/>
</p>

Mitofinder is a pipeline to assemble mitochondrial genomes and extract mitochondrial genes from trimmed 
sequencing read data.

# Requirements

This software is suitable for all linux-like systems with gcc installed (Unfortunately, not MAC and Windows < v.10).

# Table of content

1. [Installation guide for MitoFinder](#installation-guide-for-mitofinder)
2. [How to use MitoFinder](#how-to-use-mitofinder)
3. [Detailed options](#detailed-options)
4. [INPUTS](#inputs)
5. [OUTPUTS](#outputs)
6. [UCE annotation](#uce-annotation)
7. [Associated publication](#associated-publication)
8. [How to get reference mitochondrial genomes from ncbi](#how-to-get-reference-mitochondrial-genomes-from-ncbi)
9. [How to submit your mitochondrial genome(s) to GenBank NCBI](#how-to-submit-your-mitochondrial-genome(s\)-to-genbank-ncbi)

# Installation guide for MitoFinder

## Get MitoFinder

Clone mitofinder from [GitHub](https://github.com/RemiAllio/MitoFinder)

```shell 
git clone https://github.com/RemiAllio/MitoFinder.git
cd MitoFinder
./install.sh

PATH/TO/MITOFINDER/mitofinder -h  
```

or download [master.zip](https://github.com/RemiAllio/MitoFinder/archive/master.zip)  

```shell
wget https://github.com/RemiAllio/MitoFinder/archive/master.zip
unzip master.zip
mv MitoFinder-master MitoFinder
cd MitoFinder
./install.sh

PATH/TO/MITOFINDER/mitofinder -h  
```

## Add mitofinder to your path -> linux

```shell
cd PATH/TO/MITOFINDER/
p=$(pwd)
echo -e "\n#Path to mitofinder \nexport PATH=$PATH:$p" >> ~/.bashrc 
source ~/.bashrc  
```
  
WARNING: If you previously installed MitoFinder on your system and want to install a new version, you must replace the old MitoFinder PATH by the updated one in your ~/.bashrc file. To do so, you should edit your ~/.bashrc file, remove the lines that add MitoFinder to the PATH, and close your terminal. Then, you should open a new terminal and re-execute the command lines from above.

To check if the right version of MitoFinder is actually in your PATH:  
```shell
mitofinder -v
```
  
# How to use MitoFinder

### Assembler
First, you can choose the assembler using the following options:  
-- megahit 				(default: faster)  
-- metaspades			(recommended: a bit slower but more efficient (see associated paper). WARNING: Not compatible with single-end reads)  
-- idba  

## Mitochondrial genome assembly  

TIP: use mitofinder --example to print usual examples of use  

### Trimmed paired-end reads
```shell
mitofinder -j [seqid] -1 [left_reads.fastq.gz] -2 [right_reads.fastq.gz] -r [genbank_reference.gb] -o [genetic_code] -p [threads] -m [memory]   
```

### Trimmed single-end reads
```shell
mitofinder -j [seqid] -s [SE_reads.fastq.gz] -r [genbank_reference.gb] -o [genetic_code] -p [threads] -m [memory]
```

### MitoFinder can be used with your own assembly (one or several contig.s in fasta format)
```shell
mitofinder -j [seqid] -a [assembly.fasta] -r [genbank_reference.gb] -o [genetic_code] -p [threads] -m [memory]
```

### Restart
Use the same command line.  
WARNING: If you want to make the assembly again (for example because it failed) you have to remove the result assembly directory. If not, MitoFinder will skip the assembly step.  

## Test case  
```shell
cd PATH/TO/MITOFINDER/test_case/
mitofinder -j Aphaenogaster_megommata_SRR1303315 -1 Aphaenogaster_megommata_SRR1303315_R1_cleaned.fastq.gz -2 Aphaenogaster_megommata_SRR1303315_R2_cleaned.fastq.gz -r reference.gb -o 5 -p 5 -m 10   
```


# Detailed options

```shell
usage: mitofinder [-h] [--megahit] [--idba] [--metaspades] [-j PROCESSNAME]
                  [-1 PE1] [-2 PE2] [-s SE] [-a ASSEMBLY] [-m MEM]
                  [-l SHORTESTCONTIG] [-p PROCESSORSTOUSE] [-r REFSEQFILE]
                  [-e BLASTEVAL] [-n NWALK] [--ignore] [--out_gb]
                  [--blastidentitynucl BLASTIDENTITYNUCL]
                  [--blastidentityprot BLASTIDENTITYPROT]
                  [--blastsize ALIGNCUTOFF] [--circularsize CIRCULARSIZE]
                  [--circularoffset CIRCULAROFFSET] [-o ORGANISMTYPE] [-v]
                  [--example]

Mitofinder is a pipeline to assemble and annotate mitochondrial DNA from
trimmed sequencing reads.

optional arguments:
  -h, --help            show this help message and exit
  --megahit             Use Megahit for assembly. (Default)
  --idba                Use IDBA-UD for assembly.
  --metaspades          Use MetaSPAdes for assembly.
  -j PROCESSNAME, --seqid PROCESSNAME
                        Sequence ID to be used throughout the process
  -1 PE1, --Paired-end1 PE1
                        File with forward paired-end reads
  -2 PE2, --Paired-end2 PE2
                        File with reverse paired-end reads
  -s SE, --Single-end SE
                        File with single-end reads
  -a ASSEMBLY, --assembly ASSEMBLY
                        File with your own assembly
  -m MEM, --max-memory MEM
                        max memory to use in Go (MEGAHIT or MetaSPAdes)
  -l SHORTESTCONTIG, --length SHORTESTCONTIG
                        Shortest contig length to be used (MEGAHIT). Default =
                        100
  -p PROCESSORSTOUSE, --processors PROCESSORSTOUSE
                        Number of threads Mitofinder will use at most.
  -r REFSEQFILE, --refseq REFSEQFILE
                        Reference mitochondrial genome in GenBank format
                        (.gb).
  -e BLASTEVAL, --blaste BLASTEVAL
                        e-value of blast program used for contig
                        identification and annotation. Default = 0.000001
  -n NWALK, --nwalk NWALK
                        Maximum number of codon steps to be tested on each
                        size of the gene to find the start and stop codon
                        during the annotation step. Default = 200 (600 bases)
  --ignore              This option tells MitoFinder to ignore the non-
                        standart mitochondrial genes.
  --out_gb              Do not create annotation output file in GenBank
                        format.
  --blastidentitynucl BLASTIDENTITYNUCL
                        Nucleotide identity percentage for a hit to be
                        retained. Default = 50
  --blastidentityprot BLASTIDENTITYPROT
                        Amino acid identity percentage for a hit to be
                        retained. Default = 40
  --blastsize ALIGNCUTOFF
                        Percentage of overlap in blast best hit to be
                        retained. Default = 30
  --circularsize CIRCULARSIZE
                        Size to consider when checking for circularization.
                        Default = 45
  --circularoffset CIRCULAROFFSET
                        Offset from start and finish to consider when looking
                        for circularization. Default = 200
  -o ORGANISMTYPE, --organism ORGANISMTYPE
                        Organism genetic code following NCBI table (integer):
                        1. The Standard Code 2. The Vertebrate Mitochondrial
                        Code 3. The Yeast Mitochondrial Code 4. The Mold,
                        Protozoan, and Coelenterate Mitochondrial Code and the
                        Mycoplasma/Spiroplasma Code 5. The Invertebrate
                        Mitochondrial Code 6. The Ciliate, Dasycladacean and
                        Hexamita Nuclear Code 9. The Echinoderm and Flatworm
                        Mitochondrial Code 10. The Euplotid Nuclear Code 11.
                        The Bacterial, Archaeal and Plant Plastid Code 12. The
                        Alternative Yeast Nuclear Code 13. The Ascidian
                        Mitochondrial Code 14. The Alternative Flatworm
                        Mitochondrial Code 16. Chlorophycean Mitochondrial
                        Code 21. Trematode Mitochondrial Code 22. Scenedesmus
                        obliquus Mitochondrial Code 23. Thraustochytrium
                        Mitochondrial Code 24. Pterobranchia Mitochondrial
                        Code 25. Candidate Division SR1 and Gracilibacteria
                        Code
  -v, --version         Version 1.2
  --example             Print getting started examples
```

# INPUTS

Mitofinder needs several files to run depending on the method you have choosen (see above):  
- [x] **Reference_file.gb**				containing at least one mitochondrial genome of reference extracted from [NCBI](https://www.ncbi.nlm.nih.gov/)
- [ ] **left_reads.fastq.gz**				containing the left reads of paired-end sequencing  
- [ ] **right_reads.fastq.gz**				containing the right reads of paired-end sequencing  
- [ ] **SE_reads.fastq.gz** 				containing the reads of single-end sequencing  
- [ ] **assembly.fasta**				containing the assembly on which MitoFinder have to find and annotate mitochondrial contig.s   

# OUTPUTS 
### Result folder  

Mitofinder returns several files for each mitochondrial contig found:  
- [x] **[job_name]_final_genes_NT.fasta**				containing the nucleotides sequences of the final genes selected from all contigs by MitoFinder   
- [x] **[job_name]_final_genes_AA.fasta**				containing the amino acids sequences of the final genes selected from all contigs by MitoFinder   
- [x] **[job_name]_mtDNA_contig.fasta**				containing a mitochondrial contig  
- [x] **[job_name]_mtDNA_contig.gff**				containing the final annotation for a given contig (GFF3 format) 
- [x] **[job_name]_mtDNA_contig.tbl**				containing the final annotation for a given contig (Genbank submission format)
- [x] **[job_name]_mtDNA_contig.gb** 				containing the final annotation for a given contig (Genbank format for visualization)
- [x] **[job_name]_mtDNA_contig_genes_NT.fasta** 				containing the nucleotide sequences of annotated genes for a given contig    
- [x] **[job_name]_mtDNA_contig_genes_AA.fasta** 				containing the amino acids sequences of annotated genes for a given contig    
- [x] **[job_name]_mtDNA_contig.png** 				schematic representation of the annotation of the mtDNA contig    


# UCE annotation
MitoFinder starts by assembling both mitochondrial and nuclear reads. It is only in a second step that mitochondrial contigs are identified and extracted.
MitoFinder thus provides UCE contigs already assembled and the annotation can be done from the following file:  
- **[job_name]link_[assembler].scafSeq** 	containing all assembled contigs from raw reads. 

To do so, we recommend the PHYLUCE pipeline, which is specifically designed to annotate ultraconserved elements (Faircloth  2015; Tutorial: https://phyluce.readthedocs.io/en/latest/tutorial-one.html#finding-uce-loci).  
You can thus use the file **[job_name]link_[assembler].scafSeq** and start the pipeline at the **"Finding UCE"** step.  
  
# Associated publication  
  
If you use MitoFinder, please cite:  
  
Allio, R., Schomaker-Bastos, A., Romiguier, J., Prosdocimi, F., Nabholz, B., & Delsuc, F. (2019). **MitoFinder**: efficient automated large-scale extraction of mitogenomic data in target enrichment phylogenomics. BioRxiv, 685412. https://doi.org/10.1101/685412    
  
Please also cite the following studies depending on the option chosen for the assembly step in MitoFinder:    
  
Nurk, S., Meleshko, D., Korobeynikov, A., & Pevzner, P. A. (2017). **metaSPAdes**: a new versatile metagenomic assembler. Genome research, 27(5), 824-834.  
Li, D., Luo, R., Liu, C. M., Leung, C. M., Ting, H. F., Sadakane, K., ... & Lam, T. W. (2016). **MEGAHIT v1.0**: a fast and scalable metagenome assembler driven by advanced methodologies and community practices. Methods, 102, 3-11.  
Peng, Y., Leung, H. C., Yiu, S. M., & Chin, F. Y. (2012). **IDBA-UD**: a de novo assembler for single-cell and metagenomic sequencing data with highly uneven depth. Bioinformatics, 28(11), 1420-1428.  
 
For the tRNA annotation:  
  
Laslett, D., & Canbäck, B. (2008). **ARWEN**: a program to detect tRNA genes in metazoan mitochondrial nucleotide sequences. Bioinformatics, 24(2), 172-175.  
  
# HOW TO GET REFERENCE MITOCHONDRIAL GENOMES FROM NCBI  

1. Go to [NCBI](https://www.ncbi.nlm.nih.gov/)  
2. Select "Nucleotide" in the search bar  
3. Search for mitochondrion genomes:  
- [x] RefSeq (if available)  
- [x] Sequence length from 12000 to 20000  
4. Download complete record in GenBank format  

![](/image/NCBI.png)

# How to submit your mitochondrial genome(s) to GenBank NCBI   

## Submission with BankIt

If you have few mitochondrial genomes to submit, you should be able to do it with [BankIt](https://submit.ncbi.nlm.nih.gov/about/bankit/) through the NCBI submission portal.  


## Submission with tbl2asn

If you want to submit several complete or partial mitogenomes, we designed MitoFinder to facilitate the submission using [tbl2asn](https://www.ncbi.nlm.nih.gov/genbank/tbl2asn2/).  
tbl2asn requires:  
- [x] **Template file**				containing a text ASN.1 Submit-block object (suffix .sbt). [Create submission template](https://submit.ncbi.nlm.nih.gov/genbank/template/submission/).   
- [x] **Nucleotide sequence data**				containing the mitochondrial sequence(s) and associated information (suffix .fsa).  
- [x] **Feature Table**				containing annotation information for the mitochondrial sequence(s).  
- [] **Comment file**				containing assembly and annotation method information (assembly.cmt). [Create comment template](https://submit.ncbi.nlm.nih.gov/structcomment/nongenomes/)  
  
### Creating a compatible FASTA file

Because tbl2asn requires the FASTA file to contain informations associated with the data, we wrote a script to create a FASTA file containing the mitochondrial contig(s) found by MitoFinder for each species (Seq_ID) with the associated information.
This script and the associated example files can be found in the MitoFinder directory named "NCBI_submission".

#### INPUT
- [x] index_file.csv				A CSV file (comma-delimited table) containing the metadata information.

The headers of the index file are as follows:
Directory path, Seq ID, organism, location, mgcode, SRA Run1, SRA Run2, keywords ...

The first two columns are mandatory and the names cannot be changed but you can complete the index file with the different [source modifiers](https://www.ncbi.nlm.nih.gov/Sequin/modifiers.html) of NCBI by adding columns in the index file.  
NOTE: SRA Run1 = SRRXXXXXXXX1 (corresponding to SRRXXXXXXXX.1); SRA Run2 = SRRXXXXXXXX2. These columns could be left empty if your mitogenome is not associated with a particular SRA accession. However, if you use them, these two column names cannot be changed!

The directory path correponds to the path where the [Seq_ID]_mtDNA_contig.fasta file, or [Seq_ID]_mtDNA_contig_\*.fasta files if you have several contigs for the same individual, could be found. If left blank, the script will search for the contig in the directory where you run the script from (./).

TIPS:  
(1) You can copy or link (symbolic links) all your FASTA and TBL contig files in the same directory and run the script from this directory.  
(2) You can leave blanks in the index file if some species do not need a given source modifier.  
 
```shell 
/PATH/TO/MITOFINDER/NCBI_submission/create_tbl2asn_files.py -i index_file.csv
```


#### OUTPUT
- [x] **[Seq_ID].fsa**				new FASTA file containing all mtDNA contigs and the information for a given [Seq_ID]
- [x] **[Seq_ID].tbl**				new TBL file containing all mtDNA contigs and the information for a given [Seq_ID]


### Command line to run tbl2asn

Once your FASTA and TBL files have been created, you can run [tbl2asn](https://www.ncbi.nlm.nih.gov/genbank/tbl2asn2/) (download [here](ftp://ftp.ncbi.nih.gov/toolbox/ncbi_tools/converters/by_program/tbl2asn/) as follows:

```shell 
tbl2asn -t template.sbt -i [Seq_ID].fsa -V v -w assembly.cmt -a s
```

This command will create several files:
- [x] **[Seq_ID].sqn**				Submission file (.sqn) to be sent by e-mail to gb-sub@ncbi.nlm.nih.gov
- [x] **[Seq_ID].val**				Containing ERROR and WARNING values associated with tbl2asn. (ERROR explanations [here](https://www.ncbi.nlm.nih.gov/genbank/genome_validation/)

If you don't have any error and you are happy with the annotation, you can submit your mitochondrial contig(s) by sending the .sqn files to gb-sub@ncbi.nlm.nih.gov


