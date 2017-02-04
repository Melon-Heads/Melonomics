#!/usr/bin/env/ python

import sys
sys.path.append('/mnt/c/Users/Nadim/Downloads/biopython-1.68')

from Bio.Blast.Applications import NcbiblastnCommandline


###########################
#       Local BLAST       #
###########################



# This script REQUIRES:
# --> BLAST+ (downloadable from: https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
# --> A BLAST database (using nt for nucleotides), downloadable from: ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/.
# --> BioPython

# Once the database has been downloaded, it must be formatted:
#              formatdb - i nt -p F -V
# -p F [Format according to nucleotides]
# -V [Verbose]


# ------------------ 
# | BLAST Function | 
# ------------------ 


#def blastFastaFile(fastaFile):
       # blastall -i fastaFile -d nt -p blastn -b 1 -m 7 -o fastaFile"_Output.xml"


def fastaBLAST(fastaFile):
	blastn_cline = NcbiblastnCommandline(query=fastaFile, db="nt", evalue=0.001, outfmt=6, out=fastaFile+"_Output.txt", max_target_seqs=1)
	blastn_cline
	print(blastn_cline)
	stdout, stderr = blastn_cline()


# Parameters for blastall:
# --> -i or query [Input fasta file]
# --> -d or db [Database to compare against]
# --> -p [Type of blast carried out (blastn for nucleotides)]
# --> -b or max_target_seqs [Number of hits to show]
# --> -m [Format of the output file (Non-command line: 7 = XML, 8 = Tabular; Command line = 5, 6 respectively)]
# --> -o [Name of output file]


# The function is called with files passed through:
fastaBLAST('0hourR1.fasta')
fastaBLAST('24hourR1.fasta')
fastaBLAST('8hourR1.fasta')
