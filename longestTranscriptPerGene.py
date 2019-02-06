#!/usr/bin/python

import sys, getopt, pandas, argparse, warnings
import numpy as np
from gtfparse import read_gtf
warnings.filterwarnings("ignore")

parser = argparse.ArgumentParser(description="Script to identify the longest transcript isoform per gene and print it.")
parser.add_argument("input",help="Input GTF file.")
parser.add_argument("output",help="Output GTF file (will contain only longest transcript per gene.)")
args=parser.parse_args()

df= read_gtf(args.input)
f= open(args.output,'w')

#Create df of all exon lines from gtf
df_exons = df[df["feature"] == "exon"]
#Create df of all transcripts from gtf
df_transcripts = df[df["feature"] == "transcript"]


#create list of genes
df_genes = np.unique(np.array(df_transcripts.gene_id.values.tolist()))
    
#Loop over every gene
for i in df_genes:
    #Create list of transcripts for the gene
    transcripts=(df_transcripts[df_transcripts["gene_id"] == i].transcript_id.values.tolist())
    longest_trans = ""
    longest_trans_length = 0
    #loop over every transcript and work out the length
    for j in transcripts:
        curTranscript=df_exons[df_exons["transcript_id"] == j]
        transLen=(sum(curTranscript.end.values-curTranscript.start.values))
        if (transLen > longest_trans_length):
            longest_trans_length=transLen
            longest_trans = j
    #print (longest_trans , longest_trans_length)
    transcript_gtf=(df[df["transcript_id"] == longest_trans ])
   #print(transcript_gtf)
    for index, row in transcript_gtf.iterrows():
        #The following syntax doesnt work due to the need for differnt spacings between attributes in the 9th field. Therefore need the long solution below.
        #gtf1=(row['seqname'], row['source'], str(row['feature']), row['start'], row['end'], row['score'], row['strand'], row['frame'], "gene_id ")
        #print ('\t'.join(map(str,gtf)))
        f.write (row['seqname'])
        f.write ("\t")
        f.write (row['source'])
        f.write ("\t")
        f.write (row['feature'])
        f.write ("\t")
        f.write (str(row['start']))
        f.write ("\t")
        f.write (str(row['end']))
        f.write ("\t")
        f.write (str(row['score']))
        f.write ("\t")
        f.write (row['strand'])
        f.write ("\t")
        f.write (str(row['frame']))
        f.write ("\t gene_id \"")
        f.write (row['gene_id'])
        f.write ("\"; transcript_id \"")
        f.write (row['transcript_id'])
        if (row['gene_name'] != ""):
            f.write ("\"; gene_name \"")
            f.write (row['gene_name'])
        if (row['ref_gene_id'] != ""):
            f.write ("\"; ref_gene_id \"")
            f.write (row['ref_gene_id'])
        if (row['exon_number'] != ""):
            f.write ("\"; exon_number \"")
            f.write (row['exon_number'])
        f.write ("\";\n")
