#!/usr/bin/env python

import sys, re, csv, itertools # This was used just to add biopython packages to the python pat$
sys.path.append('/mnt/c/Users/Nadim/Downloads/biopython-1.68')
import pandas as pd


############################################################
# Obtaining the Gene IDs and FPKM Values from BLAST Output #
############################################################

samDict = {} # This needs to be global to add multiple blast results.

# Function to create dictionaries of the gene IDs and FPKM value per sample.
def createDict(blastOut):
	global samDict
	inputBlast = open(blastOut, "r") # Read in the file passed through the function.

	# Lists to append gene IDs and FPKM values to.
	genIDs = []
	FPKMvals = []
	
	# Dictionary of gene and FPKM values.
	genFPKM = {}

	for word in inputBlast:
		getGen = re.findall('[a-zA-Z]\w\w+\.\w+', str(word)) # Regular expression used to obtain the gene IDs. 
		getFPKM = re.findall('\d+\.\d+', str(word)) # Regular expression to obtain the FPKM values.

		genIDs.append(getGen)
		FPKMvals.append(getFPKM)

		genAcc = [l[0] for l in genIDs] # Remove any lists within the list.
		FPKM = [l[0] for l in FPKMvals] # Removes any lists within the list.

	genFPKM.update(zip(genAcc, FPKM))
	# Join the gene ID and FPKM dictionary to a separate dictionary which specifies each sample.
	samDict.update({blastOut: genFPKM}) 
	#print samDict



#####################################################################
# Creating a Sample, Gene ID and FPKM Matrix and Saving to CSV File #
#####################################################################


def dictToMatrix(dict):
	mat = pd.DataFrame(dict)
	matCSV = mat.fillna(0) # This converts empty cells to display 0.
	matCSV.to_csv('Gene_FPKM_Sample.csv')



# Passing filenames of BLAST output to the function, creating a dictionary.
# The dictionary can be used to generate a matrix.
# This matrix can be saved to a CSV file.
createDict('0hourR1.fasta_Output.txt')
createDict('8hourR1.fasta_Output.txt')
createDict('24hourR1.fasta_Output.txt')

dictToMatrix(samDict)

