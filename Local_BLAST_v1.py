#!/usr/bin/env python


##########################
#      Local BLAST       #
##########################

# Make sure you have BLAST+ downloaded from https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download
# A quick way of BLASTing and writing into an XML output file of the results
# First download the download the database - nt for nucleotides. Then using formatdb, format the database to use for the BLAST.
# The database can be downloaded from NCBI's FTP site of databases (ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/) and choose the nt.gz file to download.
# Alternatively you can use the wget command with the URL.

#  Once downloaded, the database must be formatted using the following command:
formatdb -i nt -p F -V # -p F is for nucleotides and -V means verbose.

# The function to BLAST a fasta file.
# Use one blastall function, NOT both, whichever suits your fancy :)
def blastFastaFile (fastaFile):
	blastall -i fastafile -d nt -p blastn -b 10 -m 8 > blast.out
	blastall -i Example.fasta -d nt -p blastn -b 1 -m 7 -o "Local_Blast_Output.xml" # This produces an XML file instead.
# Parameters for blastall:
# -i = The input fasta file. This can be changed to a variable.
# -d = The type of database to be compared against, described above.
# -p = The type of blast carried out - blastn is for nucleotides.
# -b = The number of results you want to show.
# -m = The format of the output file (7 = XML, 8 = Tabular)
# -o = The name of the output file.
# blast.out is the output file.


# Ignore the following:

############################
# The best way of BLASTing #
############################


#import sys
#sys.path.append('/mnt/c/Users/Nadim/Downloads/biopython-1.68')
#from Bio.Blast.Applications import NcbiblastxCommandline
# help(NcbiblastxCommandline)

#blastx_cline = NcbiblastxCommandline(query="Example.fasta", db="nt", evalue=0.001, outfmt=7, out="output.xml")
#blastx_cline
#print(blastx_cline)
#stdout, stderr = blastx_cline()
