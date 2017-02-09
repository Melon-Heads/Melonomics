import os
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
	return send_file("/mnt/c/Users/Nadim/Documents/QMUL_Level7/Group_Software_Project/Melonomics/flask/venv/HCA.pdf", attachment_filename="HCA.pdf")


####################################################
@app.route("/analyse/")
def index():
    return render_template("upload.html")

@app.route("/upload", methods=['POST'])
def upload():
    target = os.path.join(APP__ROOT, 'data/')
    print(target)

    if not os.path.isdir(target):
        os.mkdir(target)

    for file in request.files.getlist("file"):
        print(file)
        filename = file.filename
        destination = "/".join([target, filename])
        print(destination)
        file.save(destination)

    return render_template("upload2.html")

######################################################
@app.route("/jss", methods=['POST'])
def jss():
    target = os.path.join(APP__ROOT, 'data/')
    print(target)

    if not os.path.isdir(target):
        os.mkdir(target)

    for file in request.files.getlist("file"):
        print(file)
        filename = file.filename
        destination = "/".join([target, filename])
        print(destination)
        file.save(destination)
    return render_template("upload3.html")
########################################################
@app.route("/blast", methods=['POST'])
def blast():
    target = os.path.join(APP__ROOT, 'data/')
    print(target)

    if not os.path.isdir(target):
        os.mkdir(target)

    for file in request.files.getlist("file"):
        print(file)
        filename = file.filename
        destination = "/".join([target, filename])
        print(destination)
        file.save(destination)
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

