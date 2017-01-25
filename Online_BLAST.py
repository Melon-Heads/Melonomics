#!/usr/bin/env python

import sys
sys.path.append('/mnt/c/Users/Nadim/Downloads/biopython-1.68') # Required to bring Biopython into the python environment.

from Bio.Blast import NCBIWWW, NCBIXML




##########################
# BLAST using Query File #
##########################


fasta_string = open("noHeaderSample.fasta").read()
result_handle = NCBIWWW.qblast("blastn", "nt", fasta_string) 

save_file = open("myBlast_v2.xml", "w") 
save_file.write(result_handle.read())
save_file.close()

result_handle.close()
result_handle = open("myBlast_v2.xml") # To open the output file.




############################
# Parsing the BLAST Output #
############################


blast_records = NCBIXML.read(result_handle) # .read is required for a single sequence. For multiple sequences, use .parse instead.

eValThres = 0.04
for alignment in blast_record.alignments:
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
			print(hsp.match[0.75] + dots)
			print(hsp.sbjct[0:75] + dots)

