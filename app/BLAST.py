#!/usr/bin/env/ python

import sys, re, csv, os, glob

# The following is for those who need to append biopython packages to the Python PATH:
sys.path.append('/mnt/c/Users/Nadim/Downloads/biopython-1.68')

import pandas as pd
from Bio.Blast.Applications import NcbiblastnCommandline


###########################
#       Local BLAST       #
###########################


# This script REQUIRES:
# --> BLAST+ (downloadable from: https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
# --> A BLAST database (using nt for nucleotides), downloadable from: ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/.
# --> Biopython

# If a database is downloaded, it must be formatted:
#              formatdb - i nt -p F -V
# -p F [Format according to nucleotides]
# -V [Verbose]


# ------------------
# | BLAST Function |
# ------------------


def fastaBLAST(fastaFile):
        blastn_cline = NcbiblastnCommandline(query=fastaFile, db="/mnt/c/Users/Nadim/Documents/QMUL_Level7/Group_Software_Project/Scripts/nt", evalue=0.001, outfmt=6, out=fastaFile+"_Output.txt", max_target_seqs=1, num_threads=4)
        blastn_cline
        print(blastn_cline)
        stdout, stderr = blastn_cline()


# Parameters for blastall and command line blast:
# --> -i or query [Input fasta file]
# --> -d or db [Path to Database to compare against]
# --> -p [Type of blast carried out (blastn for nucleotides)]
# --> -b or max_target_seqs [Number of hits to show]
# --> -m [Format of the output file (Non-command line: 7 = XML, 8 = Tabular; Command line = 5, 6 respectively)]
# --> -o [Name of output file]
# --> num_threads [Numbers of CPUs utilised for the search]



# ------------------------------------------------------------
# | Obtaining the Gene IDs and FPKM Values from BLAST Output |
# ------------------------------------------------------------


samDict = {} # This needs to be global to add multiple blast results.
vecList = [] # A list to add the sample code for class vectors.

# Function to create dictionaries of the gene IDs and FPKM value per sample.
def createDict(blastOut, sampleCode):
        global samDict, vecList
        inputBlast = open(blastOut, "r") # Read in the file passed through the function.

        # Lists to append gene IDs and FPKM values to.
        genIDs = []
        FPKMvals = []
	# Append the sample code to the globally defined list.
	vecList.append(sampleCode)

        # Dictionary of gene and FPKM values.
        genFPKM = {}

        for word in inputBlast:
                getGen = re.findall('[a-zA-Z]\w\w+\.\w+', str(word)) # Regular expression used to obtain the gene IDs.
                genIDs.append(getGen)

		m = re.match(r"asmbl_\d+;(\d+\.\d+)", str(word)) # Regular expression to obtain the FPKM values.
		if m:
			getFPKM = m.group(1)
			FPKMvals.append(getFPKM)
				
                genAcc = [l[0] for l in genIDs] # Remove any lists within the gene ID list.
        
	# Join the gene ID and corresponding FPKM values
	genFPKM.update(zip(genAcc, FPKMvals))

        # Join the gene ID and FPKM dictionary to a separate dictionary which specifies each sample.
        samDict.update({blastOut: genFPKM})


# ---------------------------------------------------------------------
# | Creating a Sample, Gene ID and FPKM Matrix and Saving to CSV File |
# ---------------------------------------------------------------------


def dictToMatrix(dict):
	
	# A normal matrix from the dictionary formed.
        mat = pd.DataFrame(dict)
	matCSV = mat.fillna(0) # This converts empty cells to display 0.
	
	# Another matrix with the class vector mentioned in the second row.		
	vec = pd.DataFrame(dict)
        vec.loc[1] = vecList	
	vec = vec.sort()
	vecCSV = vec.fillna(0)

	# Create the CSV files.
        matCSV.to_csv('Gene_FPKM_Sample.csv')
	vecCSV.to_csv('Class_Vector.csv')


# ===================================================================== #



# To call the functions, the correct files must be chosen, found in the specific folders:
# The fasta files are sent to the BLAST function.
# The BLAST output files are utilised to create a dictionary which will be used to make a matrix.

os.chdir("/mnt/c/Users/Nadim/Documents/QMUL_Level7/Group_Software_Project/Melonomics/flask/venv/data/CTRL")
for file in glob.glob("*.fasta"):
	fastaBLAST(file)
for file in glob.glob("*Output.txt"):
	CTRLcode = 1
	createDict(file, CTRLcode)


os.chdir("/mnt/c/Users/Nadim/Documents/QMUL_Level7/Group_Software_Project/Melonomics/flask/venv/data/DS1")
for file in glob.glob("*.fasta"):
	fastaBLAST(file)
for file in glob.glob("*Output.txt"):
	DS1code = 2
	createDict(file, DS1code)


os.chdir("/mnt/c/Users/Nadim/Documents/QMUL_Level7/Group_Software_Project/Melonomics/flask/venv/data/DS2")
for file in glob.glob("*.fasta"):
        fastaBLAST(file)
for file in glob.glob("*Output.txt"):
	DS2code = 3
	createDict(file, DS2code)


dictToMatrix(samDict)


