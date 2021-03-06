#Python script to merge to results of HtSeq into a data matrix which can be used in R

#Importing modules
import sys
import re
import os

#Reading the input files
#file_number= len(sys.argv)-1
file_names=list(sys.argv)
del file_names[0]

#Start looping through the files to extract the gene names and counts and store into dictionarys
d={}
for file in file_names:
	d[file]={}
	for line in open(file):
		value= line
		match= re.match("^(?P<gene>[^\t]+)\t(?P<count>[^\t\n]+)\n$", value)
		if match:
			geneStr= str(match.group("gene"))
			geneCount= str(match.group("count"))
			d[file][geneStr]=geneCount

	
#Loop through each gene and print the counts for each sample into a file
f1=open('./htseqtemp.txt', 'w+')
for file in d:
	f1.write('\t' + file)
f1.write('\n')
for gene in d[d.keys()[0]]:
	f1.write(gene + '\t')
	for file in d:
		f1.write(str(d[file][gene]) + '\t')
	f1.write('\n') 

f2=open('./htseqcounts.txt', 'w+')

for line in open('htseqtemp.txt'):
	value = line
	match = re.match("^_",value)
	if match:
		continue
	else:
		f2.write(line)	

os.remove('htseqtemp.txt')
