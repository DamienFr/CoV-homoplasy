# SARS-CoV2 SNP Calling pipeline
Cluster nucleotide sequences using ANI (Average Nucleotid Identity) as distance estimator.

This tool was designed for the clustering of very large sequences (no limit in size) such as chromosome. It is based on an available ANI script (https://github.com/chjp/ANI). Scripts are written in shell bash and perl. Unix environment required. NCBI BLAST must be installed.

Usage: Launch_ANI.sh -f fasta_file [-i 90] [-l 90] [-s 1020] [-k] [-h]

[ ] means argument is optional


available arguments:

	-f	input file, fasta format
	-i	percentage identity cutoff for two DNA fragments to be considered from the same cluster (default 90)
	-l	minimum match length match between two fragment to be consider for clustering. percentage length of the longuest sequence (default 90)
	-s	segment length analysis window. This length must be shorter than the shorter studied sequence. Should be < 10% of the shorter sequence's length for best accuracy, but shouldn't be < 200
	-k	activate log mode. Keep the two temporary subfolders created for the analysis
	-h	display this help


The folder "scripts" contains 4 scripts that represent the backbone of the pipeline. Any script can be run alone without any arguments to get to its help instructions.

NOTE: you may have to modify the permission of the script using the chmod command i.e. "chmod +x Launch_ANI.sh"
