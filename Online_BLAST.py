#!/usr/bin/env python

import sys
sys.path.append('/mnt/c/Users/Nadim/Downloads/biopython-1.68') # Required to bring Biopython into the python environment.
# Just paste the path to the biopython-1.68 (or equivalent file) into the round brackets above.

from Bio.Blast import NCBIWWW, NCBIXML



##########################
# BLAST using Query File #
##########################


fasta_string = open("Example.fasta").read()
result_handle = NCBIWWW.qblast("blastn", "nt", fasta_string) # Carries out the BLAST search. 

save_file = open("myBlast_v2.xml", "w") # Opens a file to write the results to.
save_file.write(result_handle.read())
save_file.close()
result_handle.close() # Closing both input and output files.



############################
# Parsing the BLAST Output #
############################


blast_records = NCBIXML.read(result_handle) # .read is required for a single sequence. For multiple sequences, use .parse instead.

# To print the alignments:
eValThres = 0.04
for alignment in blast_records.alignments:
	for hsp in alignment.hsps:
		if hsp.expect < eValThres:
			print('****Alignment****')
			print('Sequence:', alignment.title)
			print('Length:', alignment.length)
			print('E-Value:', hsp.expect)
			if len(hsp.query) > 75:
				dots = '...'
			else:
				dots = ''
			print(hsp.query[0:75] + dots)
			print(hsp.match[0:75] + dots)
			print(hsp.sbjct[0:75] + dots)

result_handle.close()
