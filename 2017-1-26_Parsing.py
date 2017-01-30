#!/usr/bin/env python

import sys, re, csv
sys.path.append('/mnt/c/Users/Nadim/Downloads/biopython-1.68') # Required to bring Biopython into the python environment.
import numpy as np
from Bio import SeqIO
from Bio.Blast import NCBIXML

########################################
# Obtaining Gene IDs from BLAST Output #
########################################


def readBLASToutput(blastOutputFile):
	result_handle = open("myBlast_v2.xml") # Change this to blastOutputFile
	blast_records = NCBIXML.read(result_handle) # .read is required for a single sequence whereas .parse for multiple.

	eValThres = 0.04
	genList = [] # An empty list to introduce all gene accession numbers to.
	for alignment in blast_records.alignments:
        	for hsp in alignment.hsps:
                	if hsp.expect < eValThres: # Retrieves the top 5 gene accession numbers.
                        	#print alignment.title, "\n"
				Accession = re.findall(r'\w+\.\w+', str(alignment.title))
				genList.append(Accession) # Adds each accession number to the list.	
				genAcc = [l[0] for l in genList] # To remove the separate lists within the list.
	print genAcc
	result_handle.close()



#########################################
# Obtaining FPKM Values from Fasta File #
#########################################


def FPKMfromFasta(inputFastaFile):
	FPKMlist = []
	for seq_record in SeqIO.parse("Example.fasta", "fasta"): # Change to inputFastaFile.
		#print(seq_record.id)
		FPKM = re.findall(r'\d+\.\d*', seq_record.id)
		FPKMlist.append(FPKM)
		FPKMval = [l[0] for l in FPKMlist]
	print FPKMval



#################################
# Generating a CSV file of Data #
#################################


#keys = genAcc # Puts all the gene IDs as keys.
#values = FPKMval # Puts all the FPKM values specific to the gene IDs as values.
#matrixDict = dict(zip(keys, values))
#print matrixDict

genArray = np.column_stack((genAcc, FPKMval))
print genArray

with open("Gene_FPKM_Sample.csv", "w") as csvfile:
	#fieldnames = ["geneID", "Sample_1"]
	writer = csv.writer(csvfile)#, fieldnames=fieldnames)
	data = [#["geneID", "Sample_1"],
		[x for x in genAcc], [y for y in FPKMval]]	
	#writer.writeheader()
	for row in data:
		writer.writerow(row)

	#writer.writerows(data)
		#x = x + 1
	
	#x = 0
	#for key,value in matrixDict.items():
		#writer.writerow([key,value])
	#for x in matrixDict:
		#writer.writerow(matrixDict.keys(x), matrixDict.values(x))
		#x = x + 1

# THE DICTIONARY METHOD, ALTHOUGH SOUND HAS AN ISSUE WITH STRINGS, ETC.
# Can use numPy. First with the two lists, use npArray = np.column_stack((genAcc, FPKMval))
# Then to save to a csv file, npArray.tofile('foo.csv', sep=',')
