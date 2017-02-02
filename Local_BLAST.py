#!/usr/bin/env/python


###########################
#       Local BLAST       #
###########################


# This script REQUIRES:
# --> BLAST+ (downloadable from: https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
# --> A BLAST database (using nt for nucleotides), downloadable from: ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/.
# --> BioPython

# Once the database has been downloaded, it must be formatted:
#              formatdb - i nt -p F -V
# -p F [Format according to nucleotides]
# -V [Verbose]


# ------------------ 
# | BLAST Function | 
# ------------------ 
def blastFastaFile(fastaFile):
        blastall -i fastaFile -d nt -p blastn -b 1 -m 7 -o fastaFile"_Output.xml"

# Parameters for blastall:
# --> -i [Input fasta file]
# --> -d [Database to compare against]
# --> -p [Type of blast carried out (blastn for nucleotides)]
# --> -b [Number of hits to show]
# --> -m [Format of the output file (7 = XML, 8 = Tabular)]
# --> -o [Name of output file]



# The function is called with files passed through:
blastFastaFile('0hourR1.txt')
blastFastaFile('24hourR1.txt')
blastFastaFile('8hourR1.txt')
