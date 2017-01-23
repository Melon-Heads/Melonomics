#!/usr/bin/env python


##########################
#      Local BLAST       #
##########################

# A quick way of BLASTing and writing into an XML output file of the results
 
blastall -i p53.fasta.txt -p blastp -o Test_Result.txt 


# The best way of BLASTing

from Bio.Blast.Applications import NcbiblastxCommandline
# help(NcbiblastxCommandline)

blastx_cline = NcbiblastxCommandline(query="p53.fasta", db="nr", evalue=0.001, outfmt=5, out="output.xml")
blastx_cline
print(blastx_cline)
stdout, stderr = blastx_cline()
