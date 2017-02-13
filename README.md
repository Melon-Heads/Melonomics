# mastermelon
Group project 2017

Hello,

This is the working prototype of our software - Melonomics.

Melonomics enables users to input multiple fasta files of de novo sequences to identify possible genes and compare this expression within different samples. This is presented graphically through heatmaps and volcano plots. The analysis also displays genes identified within specific samples, presented in plots such as dendrograms.

To run the software, download the 'app' folder and run 'final.py'. This file contains all the app-routes.

The app folder should consist of static, templates and data folders. The data folder is also split into the type of samples that a user can upload - control (healthy), disease state 1 and disease state 2. HTML templates are found within the templates folder and CSS and JavaScript can be found within the static folder. Note there are a large number of files within this file structures.

To run the software, you will need:
  --> Python Flask (possibly with a virtual environment)
  --> BLAST+ to run the BLAST on command line.
  --> Biopython to partner the BLAST.
  --> A local NCBI BLAST database to compare query sequences against.
  --> Any other relevant packages (e.g. pip) or those stated within scripts.
  --> R
  --> R packages stated within the routput.R script.
  

===================================================================================

This is the file hierarchy (Apologies if this does not show in the correct format):

app
  L static
  |      L [Any files within the static folder]
  |
  L templates
  |         L [HTML Templates]
  |
  L data
  |    L [User-uploaded fasta files]
  |    L [Output from BLAST and R (when the app is run)]
  |
  L BLAST.py
  |
  L final.py
  
===================================================================================


Main folders:
  - Static --> Contains CSS and Javascript for flask framework and aesthetics of the server.
  - Templates --> All HTML templates for the web server.
  - Data --> Files for the data uploaded to the server. These are the input fasta files separated depending on the type of sample, each given a separate folder, CTRL (control), DS1 (disease state 1) and DS2 (disease state 2). Output from BLAST (e.g. matrix CSV files) and R (e.g. PDF files) are stored within this folder.


Main files and their functions:
  - final.py --> Contains all app routes to enable templates to show and imports individual software components.
  - BLAST.py --> Contains code for BLAST search and generation of output for R.
  - runner.py --> Required to run R commands within a python script. This was to enable embedding within flask framework.
  - routputs.R --> Contains code for R analysis of BLAST output.

We are working to improve the software at all times. Please feel free to leave feedback and recommendations. 

Thank you, we hope you enjoy our software!

Modupeh Betts, Andrew Knowles, Nadim Rahman and Madeleine Rhodes.
