# python_tools

	Various Python scripts I've written for analysis

# merge_htseq.py

	This tool is a script to pull together the results from multiple samples put through htseq into one text file which can be used in R for analysis

# laneCounter.py
	
	Biopython script to count reads per lane from illumina paired-end read files

# nucleotideFilter.py

	Biopython script to calculate the percentage of nucleotides per read in a pair and discard both reads in pair if nucleotide % threshold is exceeded

# gtfCDSextender.py
	
	Script utilising GTFparse (https://github.com/openvax/gtfparse) to read in a GTF file and extend the CDS regions by user specified amount. Useful for mapping of RiboSeq data to transcripts.

# aaPercent.py

	Calculate the amino acid % for CDS regions in transcript fasta file

#	codonContext.py

	Calculate the codon context (frequency of codon pairs/ frequency of amino acid pairs encoded by that codon pair) for CDS regions in transcript fasta file
