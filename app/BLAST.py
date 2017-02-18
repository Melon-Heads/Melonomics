#!/usr/bin/env/ python

#===============================================#
# File Description: Carries out BLAST and post- #
#                   BLAST processing functions. #
#                                               #
# Authors: Nadim, Modupeh, Maddy, Andrew        #
#===============================================#


# Note that all code specifying file paths may need to be altered.


################################
#      Importing Packages      #
################################

import sys, re, csv, os, glob

# If biopython packages need to be brought into the Python PATH:
sys.path.append('/mnt/c/Users/Nadim/Downloads/biopython-1.68')

import pandas as pd
from Bio import SearchIO
from Bio.Blast.Applications import NcbiblastnCommandline



###########################
#       Local BLAST       #
###########################


# This script also REQUIRES:
# --> BLAST+ (downloadable from: https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
# --> A BLAST database (using nt for nucleotides), downloadable from: ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/.
# --> Biopython

# If a BLAST database is downloaded, it must be formatted:
#              formatdb - i nt -p F -V
# -p F [Format according to nucleotides]
# -V [Verbose]


# ------------------ #
# | BLAST Function | #
# ------------------ #

# This function constructs the BLAST function on command line.
def fastaBLAST(fastaFile):
        blastn_cline = NcbiblastnCommandline(query=fastaFile, db="/mnt/c/Users/Nadim/Documents/QMUL_Level7/Group_Software_Project/Scripts/nt", evalue=0.001, outfmt=6, out=fastaFile+"_Output.txt", max_target_seqs=1, num_threads=4)
        blastn_cline
        print(blastn_cline)
        stdout, stderr = blastn_cline()


# BLAST Parameters (these may differ when shown on command line):
# --> -query [Input fasta file]
# --> -db [Path to Database to compare against]
# --> -evalue [E-Value threshold (set to 0.001)]
# --> -max_target_seqs [Number of hits to show (set to 1)]
# --> -outfmt [Format of the output file (Tabular = 6)]
# --> -out [Name of output file]
# --> -num_threads [Numbers of CPUs utilised for the search]



###################################
#      Post-BLAST Processing      #
###################################


# ------------------------------------------------------------ #
# | Obtaining the Gene IDs and FPKM Values from BLAST Output | #
# ------------------------------------------------------------ #


# Here, dictionaries are utilised to correlate first gene IDs to FPKM values.
# Then the sample name to the gene IDs and FPKM values.
# This final dictionary is fed into pandas to generate a matrix.
# The type of samples or sample code (e.g. 1 = healthy) are appended to a list.
# This list is also input into the matrix to incorporate class vectors.


samDict = {} # Dictionary for sample names (with their corresponding gene IDs and FPKM values).
vecList = [] # A list for sample codes (class vectors).

# Function to create dictionaries of the gene IDs and FPKM value per sample.
def createDict(blastOut, sampleCode):
        global samDict, vecList

        # Lists to append gene IDs and FPKM values to.
        genIDs = []
        FPKMvals = []

	# Append the sample code to the globally defined list.
	vecList.append(sampleCode)

        # Dictionary of gene and FPKM values.
        genFPKM = {}

	
	# Parse the gene IDs and FPKM values and add them to a dictionary:	
	qresults = SearchIO.parse(blastOut, 'blast-tab')
	for qresult in qresults:
		query = qresult.id
		genes = qresult.hits

		getGen = re.findall('\w+\_\d+', str(genes))
		getFPKM = re.findall('\d+\.\d+', str(query))

		genIDs.append(getGen)
		FPKMvals.append(getFPKM)

	genIDs = filter(None, genIDs) # Ensures that there are no empty entries.

        genAcc = [l[0] for l in genIDs] # Remove any lists within the gene ID list.
	FPKM = [l[0] for l in FPKMvals]	
        
	# Join the gene ID and corresponding FPKM values.
	genFPKM.update(zip(genAcc, FPKM))

        # Join the gene ID and FPKM dictionary to samDict to specify each sample.
        samDict.update({blastOut: genFPKM})



# ---------------------------------------------------------------------
# | Creating a Sample, Gene ID and FPKM Matrix and Saving to CSV File |
# ---------------------------------------------------------------------


# Here samDict is changed to a data-frame/matrix..
# CSV files are generated for R analysis.

def dictToMatrix(dict):
	global vecList

	# Another matrix with the class vector mentioned in the second row.
	# This was added as a list, which was formed in the createDict function.		
	mat = pd.DataFrame(dict)
        mat.loc[0] = vecList	
	mat = mat.sort()
	matCSV = mat.fillna(0) # This converts all empty cells to display 0.

	# Create the CSV files.
        matCSV.to_csv('/mnt/c/Users/Nadim/Documents/QMUL_Level7/Group_Software_Project/mastermelon/app/data/Gene_FPKM_Sample.csv')



# ===================================================================== #



# To call the functions, the correct files must be chosen, found in the specific folders:
# The fasta files are sent to the BLAST function.
# The BLAST output files are utilised to create a dictionary which will be used to make a matrix.

os.chdir("/mnt/c/Users/Nadim/Documents/QMUL_Level7/Group_Software_Project/mastermelon/app/data/CTRL")
for file in glob.glob("*.fasta"):
	fastaBLAST(file)
for file in glob.glob("*Output.txt"):
	CTRLcode = 1
	createDict(file, CTRLcode)


os.chdir("/mnt/c/Users/Nadim/Documents/QMUL_Level7/Group_Software_Project/mastermelon/app/data/DS1")
for file in glob.glob("*.fasta"):
	fastaBLAST(file)
for file in glob.glob("*Output.txt"):
	DS1code = 2
	createDict(file, DS1code)


os.chdir("/mnt/c/Users/Nadim/Documents/QMUL_Level7/Group_Software_Project/mastermelon/app/data/DS2")
for file in glob.glob("*.fasta"):
	fastaBLAST(file)
for file in glob.glob("*Output.txt"):
	DS2code = 3
	createDict(file, DS2code)


dictToMatrix(samDict)


