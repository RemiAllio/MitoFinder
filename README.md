MitoFinder version 1.0.2		22/03/2019
Author : Rémi ALLIO

Mitofinder is a pipeline to assemble mitochondrial genomes and extract mitochondrial genes from trimmed 
sequencing read data.

# Requirements

This software is suitable for all linux-like systems with gcc installed (Unfortunately, not MAC and Windows).

# Installation guide for MitoFinder

## Download mitofinder_vX.tar.gz at [www.gitlab.com]
```console  

$ tar -zxf mitofinder_vX.tar.gz
$ cd mitofinder_vX
$ ./install.sh

$ PATH/TO/MITOFINDER_VX/mitofinder -h  

```

## Add mitofinder to your path -> linux

```shell  

cd mitofinder_vX/
p=$(pwd)
echo -e "\n#Path to mitofinder \nexport PATH=$PATH:$p" >> ~/.bashrc 
source ~/.bashrc  

```

# How to use MitoFinder
TIP: use mitofinder --example to print usual examples of use

First, you can choose the assembler using the following options:  
-- megahit 				(default: faster)  
-- metaspades			(recommended: a bit slower but more efficient (see associated paper). WARNING: Not compatible with single-end reads)  
-- idba  

### For mitochondrial genome assembly 

## Trimmed paired-end reads
```console  

mitofinder -j [jobname] -1 [left_reads.fastq.gz] -2 [right_reads.fastq.gz] -r [genbank_reference.gb] -o [genetic_code] -p [threads] -m [memory]   

```

## Trimmed single-end reads
```console  

mitofinder -j [jobname] -s [SE_reads.fastq.gz] -r [genbank_reference.gb] -o [genetic_code] -p [threads] -m [memory]

```

## MitoFinder can be used with your own assembly (one or several contig.s in fasta format)
```console  

mitofinder -j [jobname] -a [assembly.fasta] -r [genbank_reference.gb] -o [genetic_code] -p [threads] -m [memory]

```

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
  
# To cite Mitofinder

Allio R., Schomaker-Bastos A., Romiguier J., Prosdocimi F., Nabholz B. & Delsuc F. (2019). Efficient automated extraction of mitogenomic data in target enrichment phylogenomics with MitoFinder. Molecular Ecology Resources (submitted).
