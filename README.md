# Immune_receptor_NETWORK-GENERATION
Generating a data framework for B-cell and T-cell receptor sequencing data using network analysis.

IMMUNE-NETWORK-GENERATION PROGRAM DOCUMENTATION
Program developed by Rachael Bashford-Rogers (2013)
At the Wellcome Trust Sanger Institute under Prof. Paul Kellam (rb520@cam.ac.uk)


1. Introduction
To date next-generation sequencing (NGS) of BCRs have primarily focused on classifying the IgHV, D and J recombination frequencies to understand the diversity of the BCR repertoire (1-7). However, computational assignment of V-D-J sequences to reference databases results in many incompletely assigned IgHV, D and J genes even when the germline alleles are known (8). This is most likely due to SHM masking the identity of the germline genes present in the NGS, or the existence of allelic variation relative to the reference IgH genes. Further, investigation of V-D-J gene usage frequencies utilizes only part of the BCR sequence diversity with important information about the V-D-J joining regions and mutational relationships not considered.
Here we propose that analysis of the BCR sequence relationships using the full BCR V- D-J sequence is more informative for human BCR repertoire analysis than V-D-J gene classification. We show that human BCR repertoire diversity can be interpreted through full V-D- J genotype diversity using BCR networks, previously shown to be an intuitive way for understanding B-cell repertoires in zebrafish (8). In such networks, the lowest level of organisation in a population of B-cells, namely independent B-cells, is represented by sparse networks whereas highly developed (connected) networks most likely result from clonal expansions of B-cells, arising through antigenic exposure or B-cell malignancies (8).
Here, each vertex represents a different sequence, and the number of identical BCR sequences defines the vertex size. Edges are created between vertices that differ by one nucleotide. Clusters are groups of interconnected vertices (9). The program described here calculates edges between unique sequences and determines vertex sizes, creating output files in formats that can directly be used in network analysis programs such as networkx (python) or igraph (R or python).


2. User’s guide

a) Installation

Install CD-HIT(10):
• Download current CD-HIT at http://bioinformatics.org/cd-hit • Unpack the file with “tar xvf cd-hit-XXX.tar.gz --gunzip”
• Change dir by “cd cd-hit-2006”
• Compile the programs by “make”

Download Immune-Network-Generation program (python script) and set up to user space: 
• Change directory for CD-HIT program to correct location:
line #36 cd_hit_directory = “cd-hit-v4.5.4-2011-03-07/" 
• Install the following python module dependencies:
sys, collections, os, operator, networkx


b) Usage

python Generate_networks_global.py <fasta file> <sample id> <output directory>

Arguments:

<fasta file> <sample id> <output directory>
Fasta file of sequences for network analysis.* sequence IDs should show sequence multiplicity. 
Unique sample ID.
The directory in which the output files will be created.


c) Example

> python Generate_networks_global_2.0.py TEST_FILES/Sequences_Example2.fasta Example2 OUTPUT_TEST_FILES/


