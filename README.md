MitoFinder version 1.0.1		22/03/2019
Author : Rémi ALLIO

Mitofinder is a pipeline to assemble mitochondrial genomes and extract mitochondrial genes from trimmed 
sequencing read data.

# Requirements

This software is suitable for all linux-like systems with gcc installed (Unfortunately, not MAC and Windows).

# Installation guide for MitoFinder

## Download mitofinder_vX.tar.gz at [www.github.com]

$ tar -zxf mitofinder_vX.tar.gz
$ cd mitofinder_vX
$ ./install.sh

$ PATH/TO/MITOFINDER_VX/mitofinder -h

## Add mitofinder to your path -> linux

$ cd mitofinder_vX/
$ p=$(pwd)
$ echo -e "\n#Path to mitofinder \nexport PATH=$PATH:$p" >> ~/.bashrc 
$ source ~/.bashrc

# How to use MitoFinder
TIP: use mitofinder --example to print usual examples of use

First, you can choose the assembler using the following options:
-- megahit 				(default: faster)
-- metaspades			(recommended: a bit slower but more efficient (see associated paper). WARNING: Not compatible with single-end reads)
-- idba

### For mitochondrial genome assembly 

## Trimmed paired-end reads
$ mitofinder -j [jobname] -1 [left_reads.fastq.gz] -2 [right_reads.fastq.gz] -r [genbank_reference.gb] -o [genetic_code] -p [threads] -m [memory]

## Trimmed single-end reads
$ mitofinder -j [jobname] -s [SE_reads.fastq.gz] -r [genbank_reference.gb] -o [genetic_code] -p [threads] -m [memory]

## MitoFinder can be used with your own assembly (one or several contig.s in fasta format)
$ mitofinder -j [jobname] -a [assembly.fasta] -r [genbank_reference.gb] -o [genetic_code] -p [threads] -m [memory]

### Restart
Use the same command line.
WARNING: If you want to make the assembly again (for example because it failed) you have to remove the result assembly directory. If not, MitoFinder will skip the assembly step.

# OUTPUT

## Result folder

Mitofinder returns several files for each mitochondrial contig found:
- [job_name]_partial_mito_1.fasta				containing a mitochondrial contig
- [job_name]_partial_mito_1.gff				containing the final annotation for a given contig
- [job_name]_partial_mito_1.arwen				containing the result of the tRNA annotation returned by the Arwen software.
- [job_name]_partial_mito_1.gb 				containing the final annotation for a given contig (option --out_gb)
- [job_name]_final_genes.fasta				containing the final genes selected from all contigs by MitoFinder 


## UCE annotation
MitoFinder starts by assembling both mitochondrial and nuclear reads. It is only in a second step that mitochondrial contigs are identified and extracted.
MitoFinder thus provides UCE contigs already assembled and the annotation can be done from the following file:
- [job_name]link.scafSeq 	containing all assembled contigs from raw reads. 

To do so, we recommend the PHYLUCE pipeline, which is specifically designed to annotate ultraconserved elements (Faircloth  2015; Tutorial: https://phyluce.readthedocs.io/en/latest/tutorial-one.html#finding-uce-loci).
You can thus use the file [job_name]link.scafSeq and start the pipeline at the "Finding UCE" step.

# Help
usage: mitofinder [-h] [--megahit] [--idba] [--metaspades] [-j PROCESSNAME]
                  [-1 PE1] [-2 PE2] [-s SE] [-a ASSEMBLY] [-m MEM]
                  [-l SHORTESTCONTIG] [-p PROCESSORSTOUSE] [-r REFSEQFILE]
                  [-e BLASTEVAL] [--out_gb]
                  [--blastidentitynucl BLASTIDENTITYNUCL]
                  [--blastidentityprot BLASTIDENTITYPROT]
                  [--blastsize ALIGNCUTOFF] [--circularsize CIRCULARSIZE]
                  [--circularoffset CIRCULAROFFSET] [-cove COVECUTOFF]
                  [-o ORGANISMTYPE] [-v] [--example]

Mitofinder is a pipeline to assemble and annotate mitochondrial DNA from
trimmed sequencing reads.

optional arguments:
  -h, --help            show this help message and exit
  --megahit             Use Megahit for assembly. (Default)
  --idba                Use IDBA-UD for assembly.
  --metaspades          Use MetaSPAdes for assembly.
  -j PROCESSNAME, --jobname PROCESSNAME
                        Job name to be used throughout the process
  -1 PE1, --Paired-end1 PE1
                        File with forward paired-end reads
  -2 PE2, --Paired-end2 PE2
                        File with reverse paired-end reads
  -s SE, --Single-end SE
                        File with single-end reads
  -a ASSEMBLY, --assembly ASSEMBLY
                        File with your own assembly
  -m MEM, --max-memory MEM
                        max memory to use in Go (Megahit or MetaSPAdes)
  -l SHORTESTCONTIG, --length SHORTESTCONTIG
                        Shortest contig length to be used in scaffolding.
                        Default = 100
  -p PROCESSORSTOUSE, --processors PROCESSORSTOUSE
                        Number of threads Mitofinder will use at most.
  -r REFSEQFILE, --refseq REFSEQFILE
                        Reference mitochondrial genome in GenBank format
                        (.gb).
  -e BLASTEVAL, --blaste BLASTEVAL
                        e-value of blast program used for contig
                        identification and annotation. Default = 0.000001
  --out_gb              Create annotation output file in GenBank format.
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
  -cove COVECUTOFF, --covecutoff COVECUTOFF
                        Cove cutoff for tRNAscan-SE. Default = 7
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
  -v, --version         Version 1.0.2
  --example             Print getting started examples
  
# To cite Mitofinder

Allio R., Schomaker-Bastos A., Romiguier J., Prosdocimi F., Nabholz B. & Delsuc F. (2019). Efficient automated extraction of mitogenomic data in target enrichment phylogenomics with MitoFinder. Molecular Ecology Resources (submitted).
