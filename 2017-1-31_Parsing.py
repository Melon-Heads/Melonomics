#!/usr/bin/env python

import sys, re, csv # This was used just to add biopython packages to the python path.
sys.path.append('/mnt/c/Users/Nadim/Downloads/biopython-1.68')
import pandas as pd
import numpy as np

from Bio import SeqIO
from Bio.Blast import NCBIXML


##########################
# Obtaining the Gene IDs #
##########################


result_handle = open("Local_Blast_Output.xml") # This is the input file to be parsed.

blast_records = NCBIXML.parse(result_handle)
blast_record = next(blast_records)
genList = []
for blast_record in blast_records:
	for alignment in blast_record.alignments:
		for hsp in alignment.hsps:
			if hsp.expect < 0.01:
				#print("*****Alignment*****")
				#print(alignment.title)
				Accession = re.findall('\w+\.\w+', str(alignment.title))
				genList.append(Accession)
				genAcc = [l[0] for l in genList] # To remove any lists that are within the list.
#print genAcc
genAcc = genAcc[0:5] # This was used to correlate to the FPKM values.
result_handle.close()

# Note that the parser only picks up hits.
# Therefore some gene IDs in the blast output do not show up here.



#########################################
# Obtaining FPKM Values from Fasta File #
#########################################


FPKMlist = []
for seq_record in SeqIO.parse("Example.fasta", "fasta"):
        FPKM = re.findall(r'\d+\.\d*', seq_record.id)
        FPKMlist.append(FPKM)
        FPKMval = [l[0] for l in FPKMlist]
#print FPKMval

FPKMval_2 = [1092.1, 44.2, 2.78, 1.02, 0.45]


#######################
# Creating the Matrix #
#######################


datMat = np.matrix(zip(genAcc, FPKMval, FPKMval_2))
print datMat


#####################
# Making a CSV File #
#####################


#zip(genAcc, FPKMval)
#with open("BLAST_Output.csv", "w") as csvfile:
	#fieldnames = ["geneID", "Sample_1"]
	#writer = csv.writer(csvfile, delimiter=',')
	#writer.writerows(zip(genAcc, FPKMval))



#############################################
# Appending New Column Data to the CSV File #
#############################################



# This method uses Pandas:
#df = pd.read_csv('BLAST_Output.csv', delimiter='\t')
#new_column = df['Sample_2'] + 1
#df['NextColumn'] = new_column
#df.to_csv('BLAST_Output_2.csv', sep='\t')

# A way that kind of works, but replicates on each column as opposed to just printing additional on each row.
#with open("BLAST_Output.csv", "r") as csvIn, open("BLAST_Output_2.csv", "w") as csvOut:
	#reader = csv.reader(csvIn, lineterminator='\n') 
	#writer = csv.writer(csvOut, lineterminator='\n', delimiter=',')
	#for row, val in zip(reader, FPKMval_2):
		#writer.writerow(row + FPKMval_2)

#raw_data = {'geneID': genAcc,
	#'Sample_1': FPKMval}
#df = pd.DataFrame(raw_data, columns = ['geneID', 'Sample_1'])
#df

#df.to_csv('Tester.csv')


#colDf = pd.read_csv('Tester.csv', delimiter='\t')
#new_column = colDf['Sample_1'] + 1
#colDf['NewColumn'] = new_column
#colDf.to_csv('Tester.csv', sep='\t')
