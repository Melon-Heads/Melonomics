######## Melonomics Software Development Project############
# for: QMUL MSc Bioinformatics
#Name: Melon Group (Modupeh Betts, Nadim Rahman, Maddy Rhodes, Andrew Knowles)
#Date: Submitted for assement 17/2/17



#imports python os package to enable access to the sub-modules such as os.path specifying file saving path
import os
#imports modules of flask that are required
from flask import Flask, render_template, request, url_for, redirect

app = Flask(__name__)

APP__ROOT = os.path.dirname(os.path.abspath(__file__))


#this returns the home page
@app.route('/')
def home():
	return render_template ('INDEX.html')

                                                    
#route for about us page
@app.route('/about/')
def about():
	return render_template ('ABOUTUS.html')

# this route returns the about software page
@app.route('/software/')
def software():
	return render_template ('software.html')

# this route handles 404 errors and returns the error.html page
@app.errorhandler(404)
def error(e):
	return render_template ('error.html')	
#This route is call upon when all 3 type files are submitted
#route calls on the BLAST pyhton script and runs analysis
@app.route("/Blast/")
def Blast(name=None):
	import BLAST
	print("Completed BLAST, waiting for Analysis!")
	return render_template("Blast.html")

# This shows a template which consists of a button to allow users to view R analysis.
@app.route("/R_Downloads/")
def R_Downloads():
	return render_template("Rdownloads.html")

# This displays a pdf of R analysis on the screen.
@app.route("/return_file/")
def returnFile():
	return send_file("/data/HCA.pdf", attachment_filename="HCA.pdf")

# This displays a pdf of the heatmap on screen.
@app.route("/return_heatmap/")
def returnHeat():
	return send_file("/data/heatmap.pdf", attachment_filename="Heatmap.pdf")

#this route ask for input files control files and save in CTRL folder
@app.route("/analyse/")
def index():
    return render_template("upload.html")

@app.route("/upload", methods=['POST'])
def upload():
    #creates the CTRL folder in this path saves file
    target = os.path.join(APP__ROOT, 'data/CTRL')
    print(target)

    if not os.path.isdir(target):
        os.mkdir(target)

    for file in request.files.getlist("file"):
        print(file)
        filename = file.filename
        destination = "/".join([target, filename])
        print(destination)
        file.save(destination)
        #after file saving returns to form for second DS1 uploads
    return render_template("upload2.html")

#This route ask for input files for 8hrs infection files and saves in DS1 folder

@app.route("/jss", methods=['POST'])
def jss():
    #creates the DS1 folder in this path and saves files
    target = os.path.join(APP__ROOT, 'data/DS1')
    print(target)

    if not os.path.isdir(target):
        os.mkdir(target)

    for file in request.files.getlist("file"):
        print(file)
        filename = file.filename
        destination = "/".join([target, filename])
        print(destination)
        file.save(destination)
        #files will be saved and the form for DS2 will be returned
    return render_template("upload3.html")


#This route ask for input files for 8hrs infection files and saves in DS2 folder in POST

@app.route("/blast", methods=['POST'])
def blast():
    #create DS2 folder if not exist and save files
    target = os.path.join(APP__ROOT, 'data/DS2')
    print(target)

    if not os.path.isdir(target):
        os.mkdir(target)

    for file in request.files.getlist("file"):
        print(file)
        filename = file.filename
        destination = "/".join([target, filename])
        print(destination)
        file.save(destination)
        #The uploaded files are then returned to the function with blast file

    return "This will call on function where the analysis takes place"

#This is an alternative method of saving files that we may try later
'''@app.route('/submit/', methods=['GET', 'POST'])
def submit():
	if request.method == 'POST':
		fasta = request.form['fasta']

		return render_template('fasta.html',title='Submit', fasta=fasta)

	return render_template('form.html')


@app.route('/upload', methods=['GET', 'POST'])
def upload():
	if request.method =='POST' and text in request.files:
		filename=text.save(request.files['text'])
		return filename
	return render_template ('uploads.html')'''



							
if __name__ == '__main__':
	app.run (debug=True)

