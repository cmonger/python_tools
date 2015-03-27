#Python script to merge to results of HtSeq into a data matrix which can be used in R

#Importing modules
import sys
import re


#Reading the input files
file_number= len(sys.argv)-1
file_names=list(sys.argv)
del file_names[0]

#Create a dictionary for each file

for file in file_names:
	file={}

#Start looping through the files to extract the gene names and counts and store into dictionarys

for file in file_names:
	d = {}
	for line in open(file):
		value= line
		match= re.match("^(?P<gene>[^\t]+)\t(?P<count>[^\t\n]+)\n$", value)
		if match:
			d[match.group("gene")] = match.group("count")
	file = d	

				
#Loop through each gene and print the counts for each sample into a file

	
	
