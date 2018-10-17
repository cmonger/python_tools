#!/usr/bin/python
import sys
from Bio import SeqIO
from itertools import izip_longest
from collections import Counter

#Usage nucleotideFilter.py {interleaved.fq}  {nucleotide%threshold}

file1 = sys.argv[1]
#Read in the interleaved fastq
readiter = SeqIO.parse(file1, "fastq")
lane1=0
lane2=0
lane3=0
lane4=0
lane5=0
lane6=0
lane7=0
lane8=0
#Parse through the file
for rec1, rec2 in izip_longest(readiter, readiter):
	#Read the lane ID
	if(int(rec1.id.split(":")[3]) == 1):
		lane1 = lane1+1
        if(int(rec1.id.split(":")[3]) == 2):
                lane1 = lane1+1
        if(int(rec1.id.split(":")[3]) == 3):
                lane1 = lane1+1
        if(int(rec1.id.split(":")[3]) == 4):
                lane1 = lane1+1
        if(int(rec1.id.split(":")[3]) == 5):
                lane1 = lane1+1
        if(int(rec1.id.split(":")[3]) == 6):
                lane1 = lane1+1
        if(int(rec1.id.split(":")[3]) == 7):
                lane1 = lane1+1
        if(int(rec1.id.split(":")[3]) == 8):
                lane1 = lane1+1

print "Lane 1 read-pairs:	" , lane1
print "Lane 2 read-pairs:	" , lane2
print "Lane 3 read-pairs:	" , lane3
print "Lane 4 read-pairs:	" , lane4
print "Lane 5 read-pairs:	" , lane5
print "Lane 6 read-pairs:	" , lane6
print "Lane 7 read-pairs:	" , lane7
print "Lane 8 read-pairs:	" , lane8

