#!/usr/bin/env python


##########################
#      Local BLAST       #
##########################

# A quick way of BLASTing and writing into an XML output file of the results
# First download the download the database - nt for nucleotides. Then using formatdb, format the database to use for the BLAST.
# The database can be downloaded from NCBI's FTP site of databases or using the wget command and pasting the URL.

#formatdb -i nt -p F -V # -p F is for nucleotides and -V means verbose.

# The function to BLAST a fasta file.
def blastFastaFile (fastaFile):
	blastall -i fastafile -d nt -p blastn -b 10 -m 8 > blast.out
	blastall -i Example.fasta -d nt -p blastn -b 10 -m 7 -o "Local_Blast_Output.xml" # This produces an XML file instead.




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
