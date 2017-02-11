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
	
#this route ask for input files control files and save in CTRL folder
@app.route("/analyse/")
def index():
    return render_template("upload.html")

@app.route("/upload", methods=['POST'])
def upload():
    #creates the CTRL folder in this path saves file
    target = os.path.join(APP__ROOT, 'data/CTRL/')
    print(target)

    if not os.path.isdir(target):
        os.mkdir(target)

    for file in request.files.getlist("file"):
        print(file)
        filename =file.filename
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

#This route ask for input files for 8hrs infection files and saves in DS2 folder
@app.route("/blast", methods=['POST'])
def blast():
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

