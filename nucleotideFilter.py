#!/usr/bin/python
import sys 
from Bio import SeqIO
from itertools import izip_longest
from collections import Counter

#Usage nucleotideFilter.py {interleaved.fq}  {nucleotide%threshold}

file1 = sys.argv[1]
Nthreshold = sys.argv[2]
outfile1= (file1.split(".")[0] + "_1_filtered.fq")
outputfile1=open(outfile1,"w")
outfile2= (file1.split(".")[0] + "_2_filtered.fq")
outputfile2=open(outfile2,"w")

#Read in the interleaved fastq
readiter = SeqIO.parse(file1, "fastq")

#Parse through the file
for rec1, rec2 in izip_longest(readiter, readiter):
	#Count occurance of each nucleotide
	counts1=Counter(rec1.seq)
	counts2=Counter(rec2.seq)
	#Count the sequence length
	seqLen1=len(rec1.seq)
	seqLen2=len(rec2.seq)
	#make empty dictionarys
	d1 = {}
	d2 = {}
	#Loop through the counts of each nucleotide and work out as a percentage of the sequence length
	for k,v in  counts1.items():
		d1[k]=  float(100/float(seqLen1)*v)
	for k,v in counts2.items():
		d2[k] = float(100/float(seqLen2)*v)
	#Check all the values are  below the threshold
        if (all (x <  Nthreshold for x in d1.values())) and (all (x <  Nthreshold for x in d2.values())):
		SeqIO.write(rec1, outputfile1, "fastq")
		SeqIO.write(rec2, outputfile2, "fastq")
outputfile1.close()
outputfile2.close()
