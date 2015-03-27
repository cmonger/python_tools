#Python script to merge to results of HtSeq into a data matrix which can be used in R

#Importing modules
import sys
import re


#Reading the input files
file_number= len(sys.argv)-1
file_names=list(sys.argv)
del file_names[0]

#
for file in file_names:
	print file 
