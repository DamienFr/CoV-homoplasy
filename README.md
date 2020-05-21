# SARS-CoV-2 SNP Calling pipeline

This repository contains scripts and command line used for the publication
"Emergence of genomic diversity and recurrent mutations in SARS-CoV-2"
L. van Dorp, M. Acman*, D. Richard*, L. P. Shaw*, C. E. Ford, L. Ormond, C. J. Owen, J. Pangae, C. C. S. Tan, F. A.T. Boshiere, A. Torres Ortiz, F. Balloux
\* equal contribution
https://doi.org/10.1016/j.meegid.2020.104351

This pipeline takes as input the raw sequencing reads (fastq) of 348 SARS-CoV-2 isolates and outputs a set of SNP and nucleotides alignements.

In the context of the publication, phylogenies are produced and homoplasies are identified.

The script 'Main.sh" is to be run in a unix command line step by step to ensure proper functioning.
You are likely to have to make some modifications to it in order to modify paths etc.

## Inputs (not provided)
- A fasta reference genome to be put in a folder named 00.ref
- Fastq files containing raw sequencing reads to be put in a folder named 01.sra_fastq

## Outputs (not provided)
- 07.filtrated_WGS_all_sites.fasta Multifasta aligned file 
- 07.filtrated_WGS_only_variant.fasta Multifasta aligned file 

## Dependencies
picard.jar MarkDuplicates (> v.2.22)
bwa mem 
samtools (> v.1.7 )
Perl core
Perl List::Util module (usually included in Perl core)
Perl File::Copy module (usually included in Perl core)

XML::LibXML
Getopt::Long

freebayes


This Readme file will be edited and completed more thoroughly in the days to come.

