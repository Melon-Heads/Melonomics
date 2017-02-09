# mastermelon
Group project 2017

Hello,

This is the working prototype of our software - Melonomics.

To run the software, download the 'app' folder and run 'final.py'. This file contains all the app-routes.

The app folder should consist of static, templates and data folders. The data folder is also split into the type of samples that a user can upload - control (healthy), disease state 1 and disease state 2. HTML templates are found within the templates folder and CSS and JavaScript can be found within the static folder. Note there are a large number of files within this file structures.

To run the software, you will need:
  --> BLAST+ to run the BLAST on command line.
  --> Biopython to partner the BLAST.
  --> A local NCBI BLAST database to compare query sequences against.

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

We are working to improve the software at all times. Please feel free to leave feedback and recommendations. 

Thank you, we hope you enjoy our software!

Modupeh Betts, Andrew Knowles, Nadim Rahman and Madeleine Rhodes.
