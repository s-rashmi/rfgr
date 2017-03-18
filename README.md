# rfgr
Repeat finder for whole genome and NGS Reads

RFGR has the following requirements and options for its use:

RFGR has the following input options for the user

**(A) For repetitive elements in genomes**
  * genome – to provide complete genome as input
  * contig – to provide assembled genome as input
  * k – to input the desired k-mer length (default: 16)
  * gff – to provide the GFF file to find nearest genes for repeats that are within 500bp of each other (optional)

**(B) For repetitive elements in Next-Generation Sequencing (NGS) reads**
  * reads – to provide single-end NGS reads as input
  * pe – to invoke paired-end settings
  * r1 and r2 – to provide the forward and reverse reads dataset respectively for paired-end sequencing
  * k – to input the desired k-mer length (recommended: 75% of read length)
  
For example,

`perl rfgr.pl -genome ecoli.txt -k 16 -gff ecoli_genes.gff`

**The output will contain**

(A) For whole genomes (complete and assembled), each line of the results file contains details of one pair of repeat locations. 

The output of RFGR for whole genomes contains the following fields:
1. Genome ID: The ID of the genome sequence from which the k-mers are generated
2. Position of origin: The position in the genome where the k-mer being reported was generated. Thus, it is also the position of one copy of the repeat. The postion is 1-offset and reported as a pair of coordinates delimiting its location in the input sequence.
3. Length of the repeat: The length of the maximal repeat at that position.
4. Reference sequence: The sequence in which the copy of the repeat is seen. Some organisms have more than one chromosome and copies of repeats are seen translocated between chromosomes.
5. Repeat pair positions: This column shows the position coordinates in the reference sequence, where the copy of the repeat is found. Positions on forward strand have ascending co-ordinates and positions on the reverse strand have descending co-ordinates.
6. Type of repeat: This column shows the type of repeat the pair of positions represents. This column is only present for complete genomes.
7. Closest gene: This column provides the gene that is closest to the pair of repeat positions, if the positions are less than 500 bases apart. From this column users can identify the gene that this pair of repeat may have a influence on gene regulation. This column is only present for complete genomes.

(B) For NGS Reads dataset, a new set of reads are produced from the original reads, where redundant repetitive reads are removed.
