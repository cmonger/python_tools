import sys
from Bio import SeqIO

def split_bylen(item, maxlen):
#Split a string into a list of strings of length(maxlen)
    return [item[ind:ind+maxlen] for ind in range(0, len(item), maxlen)]

def increaseDictCounts(key,value, dict):
#Check if an item is in a dictionary. If it is, increasment the count by given value. If not then add it and set the count to the value.
	if (str(key) in dict):
		dict[str(key)] += value
	else :
		dict[str(key)] = value

def translate(codonList):
  return str(map[codonList[0]] + map[codonList[1]])


map = {"TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
   "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
   "TAT":"Y", "TAC":"Y", "TAA":"*", "TAG":"*",
   "TGT":"C", "TGC":"C", "TGA":"*", "TGG":"W",
   "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
   "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
   "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
   "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
   "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
   "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
   "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
   "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
   "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
   "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
   "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
   "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G",}


dCodonPairCounts ={}
dAaPairCounts ={}

#loop through the records
for record in SeqIO.parse(sys.argv[1] ,"fasta"):
	if record.seq.count('N')==0:
		for i in (range(0,len(record.seq)/3)):
			codon=i*3
			#Check it isnt a stop codon in the first frame
			if (str(record.seq[codon:codon+3]) != "TAA" or "TAG" or "TGA"):
				#Sum the codon pairs in the dictionary
				codonPair= str(record.seq[codon:codon+6])
				if (len(codonPair) == 6):
					increaseDictCounts(codonPair,1,dCodonPairCounts)

#Translate the codonPair counts to aminoAcid pair counts
for codonPair,codonPairCount in dCodonPairCounts.items():
	codonList=split_bylen(codonPair,3)
	#Translate using the map dictionary
	aaPair= translate(codonList)
#	aaPair= str(map[codons[0]] + map[codons[1]])
	#Sum counts to the aa dictionary
	increaseDictCounts(aaPair,codonPairCount, dAaPairCounts)

#print dAaPairCounts

#Loop through every combination of Aas
#for codon1, aa1 in map.items():
#	for codon2,aa2 in map.items():
#		print aa1+aa2


#I have the  counts of every aa and codon pair. Now I need to find the codon context

for codonPair,codonPairCount in dCodonPairCounts.items():
	codonList=split_bylen(codonPair,3)
	aaPair= translate(codonList)
#	print codonPair, codonPairCount,aaPair, dAaPairCounts[aaPair] , float(codonPairCount)/float(dAaPairCounts[aaPair])
	print codonPair,float(codonPairCount)/float(dAaPairCounts[aaPair])

