#!/usr/bin/python

import sys, getopt, pandas, argparse, warnings
from gtfparse import read_gtf
warnings.filterwarnings("ignore")

#read sys args
parser = argparse.ArgumentParser(description="Script to modify the CDS region of transcripts in a GTF file by a specified number of nucleotides. The strandedness of the sequences is considered.")
parser.add_argument("input",help="input GTF file")
parser.add_argument("index",help="input fa index file (.fai from bowtie/samtools index)")
parser.add_argument("output",help="output GTF file (will contain only CDS regions)")
parser.add_argument("-s",'--start',default=13,type=int,help="How many nucleotides should the 5' end of the CDS be modified by?\tDefault=13")
parser.add_argument("-e",'--end',default=10,type=int,help="How many nucletides should the 3' end of the CDS be modified by?\tDefault=10")
parser.add_argument("-E",'--errorPrefix',default="gtfextender",type=str,help="File to log error output for gtf lines extending past boundaries etc.")
parser.add_argument("-c",'--cdsError',action='store_true',help="Output a list of transcripts with no CDS to a log file.")
args=parser.parse_args()

print ("\n\nExtending the CDS regions in " + args.input + " by " + str(args.start) + "nt in the start and " + str(args.end) + "nt in the end of the CDS of every transcript\n")

#Open the input and output files
df=read_gtf(args.input)
f = open(args.output, 'w')
e = open(args.errorPrefix+"CDSerror.log", 'w')
if(args.cdsError):
	c = open(args.errorPrefix+"NoCDS.log", 'w')
#Parse the gtf file to get a each of the transcripts
df_transcripts=df[df["feature"] == "transcript"]
#Create a list of all the transcript IDs
df_transcripts.list = df_transcripts.transcript_id.values.tolist()

#Loop through the fasta index to get the chr name and length
def create_index_table(index):
	df_index= pandas.read_table(index, delim_whitespace=True,header=-1)
	df_index.columns=["chromosome","length","totlength","ind1","ind2"]
	df_index.drop('totlength', axis=1,inplace=True)
	df_index.drop('ind1', axis=1, inplace=True)
	df_index.drop('ind2', axis=1, inplace=True)
	return (df_index)

def check_CDS_boundaries(transcript, chr):
	global errorCount
#	print "Transcript "+transcript +" is located on "+ chr
	df_chrlen=(df_index[df_index["chromosome"] == chr])
#	print "Chromosome " + chr + " has length " + str(int(df_chrlen.length))
	if (CDS.iloc[0,6] == "+"):
#		print "Transcript "+ transcript+ "is on the positive strand and has start and end positions " + str(startCDS)  +" "+ str(endCDS)
		if (startCDS <= 0):
			errorCount+=1
			e.write ("Dropped transcript " + transcript + " as its extended positions exceed chr boundaries (" + str(startCDS) + ' < 1)\n')
			return(1)
		elif (endCDS > int(df_chrlen.length)):
			errorCount+=1
			e.write ("Dropped transcript " + transcript + " as its extended positions exceed chr boundaries ("  + str(endCDS) + " > " + str(int(df_chrlen.length)) + ")\n")
			return(1)
		else:
			return(2)
	elif (CDS.iloc[0,6] == "-"):
#		print "Transcript "+ transcript+ "is on the negative strand and has start and end positions " + str(startCDS)  +" "+ str(endCDS)
		if (endCDS <= 0):
			errorCount+=1
			e.write ("Dropped transcript " + transcript + " as its extended positions exceed chr boundaries ("  + str(endCDS) + " < 1)\n")
			return(1)
		elif (startCDS > int(df_chrlen.length)):
			errorCount+=1
			e.write ("Dropped transcript " +transcript + " as its extended positions exceed chr boundaries (" + str(startCDS) + " > " + str(int(df_chrlen.length)) + ")\n")
			return(1)
		else:
			return(2)
	else:
			e.write ("Dropped transcript " + transcript + " as it has no strand information in the GTF input\n")
			return(1)

#A function to output the GTF file once the cds has been edited later
def return_CDS_gtf(last_CDS_row_index):
	for i in range(0,last_CDS_row_index+1):
		f.write(str(CDS.iloc[i,0]) +"\t"+ str(CDS.iloc[i,1])  +"\t"+ str(CDS.iloc[i,2]) +"\t"+ str(CDS.iloc[i,3]) +"\t"+ str(CDS.iloc[i,4]) +"\t.\t" + str(CDS.iloc[i,6]) + "\t.\tgene_id \"" + str(CDS.iloc[i,8]) +"\"; transcript_id \"" + str(CDS.iloc[i,12])+ "\"; ")
		if (str(CDS.iloc[i].gene_name) != ""):
			f.write ("gene_name \""+ str(CDS.iloc[i].gene_name + "\"; "))
		if (str(CDS.iloc[i].exon_number) != ""):
                        f.write ("exon_number \""+ str(CDS.iloc[i].exon_number + "\"; "))
		if (str(CDS.iloc[i].gene_biotype) != ""):
                        f.write ("gene_biotype \""+ str(CDS.iloc[i].gene_biotype + "\"; "))
		if (str(CDS.iloc[i].transcript_name) != ""):
                        f.write ("transcript_name \""+ str(CDS.iloc[i].transcript_name + "\"; "))
		if (str(CDS.iloc[i].transcript_biotype) != ""):
                        f.write ("transcript_biotype \""+ str(CDS.iloc[i].transcript_biotype + "\"; "))
		if (str(CDS.iloc[i].protein_id) != ""):
                        f.write ("protein_id \""+ str(CDS.iloc[i].protein_id + "\"; "))
		f.write("\n")

#Run the fasta indexer
df_index= create_index_table(args.index)

CDScount= 0
errorCount= 0
transcriptCount= 0
#Loop through each transcript
for i in range (0,(len(df_transcripts.list))):
	transcriptCount+=1
	transcript=df_transcripts.list[i]
	#Get a slice of each transcript
	df_transcript_gtf=(df[df["transcript_id"] == transcript])
	#Error check here to see if there is a CDS as some transcripts are non coding or unannotated CDS
	if ('CDS' in df_transcript_gtf.feature.values):
		CDScount=CDScount+1
		#Create new dataframe with the CDS region
		CDS=(df_transcript_gtf[df_transcript_gtf["feature"] == "CDS"])
		#If/elsif to check strandedness of the transcript annotation
		if (CDS.iloc[0,6] == "+"):
			## Modify the start location
			startCDS=CDS.iloc[0,3]-13
			CDS.iloc[0,3]=startCDS
			##Modify the end location
			last_CDS_row_index= CDS.shape[0]-1
			endCDS=CDS.iloc[last_CDS_row_index,4]+10
 			CDS.iloc[last_CDS_row_index,4]=endCDS

		elif (CDS.iloc[0,6] == "-"):
			## Modify the start location
			startCDS=CDS.iloc[0,4]+13
			CDS.iloc[0,4]=startCDS

			##Modify the end location
			last_CDS_row_index= CDS.shape[0]-1
			endCDS=CDS.iloc[last_CDS_row_index,3]-10
			CDS.iloc[last_CDS_row_index,3]=endCDS
		result = check_CDS_boundaries(transcript, df_transcript_gtf.iloc[1,].seqname)

		if (result != 1):
			return_CDS_gtf(last_CDS_row_index)
	else:
		if  (args.cdsError):
			c.write("Transcript "+ transcript + " was dropped as it has no CDS in the GTF file\n")
print  ("\n\n" + str(transcriptCount) + " transcripts were in the input GTF. " + str(CDScount) + " of these had a CDS of which " +str(CDScount - errorCount) +" have been extended successfully and output to the file " + str(args.output))
