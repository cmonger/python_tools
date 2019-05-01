#!/usr/bin/python
import sys
from Bio import SeqIO
from Bio.Alphabet.IUPAC import extended_protein
d = dict.fromkeys(list(extended_protein.letters), 0)
d["."]=0

aaCount=0

for record in SeqIO.parse(sys.argv[1], "fasta"):
	#The sequence is in the `seq` attribute of the record
#	print record.seq
	for aa in record.seq:
		d[aa] += 1
		if (aa != '.'):
			aaCount += 1

del d["."]
for key,val in d.items():
	print key, val, str((float(100)/float(aaCount)* val)) +"%"

