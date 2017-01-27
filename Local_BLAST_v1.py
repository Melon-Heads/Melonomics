#!/usr/bin/env python


##########################
#      Local BLAST       #
##########################

# A quick way of BLASTing and writing into an XML output file of the results
 

blastall -i p53.fasta.txt -d human.1.protein.faa -p blastp -o "Result.txt" # This produces the output (in a text file) in exactly the same way as online BLAST.

blastall -i p53.fasta.txt -d human.1.protein.faa -p blastp -o "Results" -m 7

blastall -i Sample.fasta -d human.1.protein.faa -p blastp -b 10 -m 8 > blast.out # Produces a blast.out file with the scores, e-values etc.
# This can be used to generate a table to show on the web page.




# The best way of BLASTing

from Bio.Blast.Applications import NcbiblastxCommandline
# help(NcbiblastxCommandline)

blastx_cline = NcbiblastxCommandline(query="p53.fasta", db="nr", evalue=0.001, outfmt=5, out="output.xml")
blastx_cline
print(blastx_cline)
stdout, stderr = blastx_cline()
