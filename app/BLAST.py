#!/usr/bin/env/ python

import sys, re, csv, os

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
# --> -d or db [Database to compare against]
# --> -p [Type of blast carried out (blastn for nucleotides)]
# --> -b or max_target_seqs [Number of hits to show]
# --> -m [Format of the output file (Non-command line: 7 = XML, 8 = Tabular; Command line = 5, 6 respectively)]
# --> -o [Name of output file]
# --> num_threads [Numbers of CPUs utilised for the search]



# ------------------------------------------------------------
# | Obtaining the Gene IDs and FPKM Values from BLAST Output |
# ------------------------------------------------------------


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



# ---------------------------------------------------------------------
# | Creating a Sample, Gene ID and FPKM Matrix and Saving to CSV File |
# ---------------------------------------------------------------------


def dictToMatrix(dict):
        mat = pd.DataFrame(dict)
        matCSV = mat.fillna(0) # This converts empty cells to display 0.
        matCSV.to_csv('Gene_FPKM_Sample.csv')



# Look inside the data folder to BLAST all contents and create a dictionary.
for file in os.listdir("/mnt/c/Users/Nadim/Documents/QMUL_Level7/Group_Software_Project/Melonomics/flask/venv/data"):
	if file.endswith(".fasta"):
		fastaBLAST(file)
		createDict(file)


# Functions to BLAST uploaded files:
#fastaBLAST('0hourR1.fasta')
#fastaBLAST('24hourR1.fasta')
#fastaBLAST('8hourR1.fasta')

# Functions to process the BLAST output:
#createDict('0hourR1.fasta_Output.txt')
#createDict('8hourR1.fasta_Output.txt')
#createDict('24hourR1.fasta_Output.txt')

dictToMatrix(samDict)


