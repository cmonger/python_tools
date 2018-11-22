
#!/usr/bin/python

import sys, getopt, pandas, argparse, warnings
from gtfparse import read_gtf
warnings.filterwarnings("ignore")

#read sys args
parser = argparse.ArgumentParser(description="Script to modify the CDS region of transcripts in a GTF file by a specified number of nucleotides. The strandedness of the sequences is considered.")
parser.add_argument("input",help="input GTF file")
parser.add_argument("output",help="output GTF file (will contain only CDS regions)")
parser.add_argument("-s",'--start',default=13,type=int,help="How many nucleotides should the 5' end of the CDS be modified by?\tDefault=13")
parser.add_argument("-e",'--end',default=10,type=int,help="How many nucletides should the 3' end of the CDS be modified by?\tDefault=10")
args=parser.parse_args()


print ("\n\nExtending the CDS regions in " + args.input + " by " + str(args.start) + "nt in the start and " + str(args.end) + "nt in the end of the CDS of every transcript\n")


#Read the input GTF file
df=read_gtf(args.input)
f = open(args.output, 'w')
#Parse the gtf file to get a each of the transcripts
df_transcripts=df[df["feature"] == "transcript"]

#A function to output the GTF file once the cds has been edited later
def return_CDS_gtf(last_CDS_row_index):
	for i in range(0,last_CDS_row_index):
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


CDScount= 0
#Loop through each transcript
for i in range(len((list(df_transcripts.transcript_id)))):
	#Get the transcript ID using the index above. 
	transcript=(df_transcripts.iloc[i,12])
	#Get a slice of the gtf file dataframe containing only lines for the current transcript in the loop
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
		return_CDS_gtf(last_CDS_row_index)

print  ("\n\n" + str(i) + " transcripts were in the input GTF. " + str(CDScount) + " of these had a CDS which has been extended and output to the file " + str(args.output))
