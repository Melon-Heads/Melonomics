import os
from flask import Flask, render_template, request, url_for, redirect
from flask.ext.uploads import UploadSet, configure_uploads, TEXT
from werkzeug.utils import secure_filename

#APP__ROOT = os.path.dirname(os.path.abspath(__file__))


app = Flask(__name__)
#give variable to upload
text = UploadSet ('text', TEXT)
#from app import views
app.config['UPLOADED_TEXT_DEST'] = 'static/uploads'
configure_uploads (app, text)

@app.route('/')
def home():
	return render_template ('index2.html')

#def index():
 #   user = {'nickname': 'Melon Heads'}  # fake user
 #   return render_template('index.html',title='Home',user=user)
                                                     
@app.route('/index/')
def guess():
	return render_template ('boot.html')

@app.route('/about/')
def about():
	return render_template ('aboutUs.html')


@app.route('/submit/', methods=['GET', 'POST'])
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
	return render_template ('uploads.html')


#########################################################
	
@app.route('/uploader/', methods = ['GET', 'POST'])
def upload_file():
   if request.method == 'POST':
      f = request.files['file']
      f.save(secure_filename(f.filename))
      return 'file uploaded successfully'
#########################################################


#@app.route("/hi")
#def index():
 #   return render_template("uploads.html")

#@app.route("/upload", methods=['POST'])
#def upload():
 #   target = os.path.join(APP__ROOT, 'data/')
  #  print(target)

   # if not os.path.isdir(target):
    #    os.mkdir(target)

    #for file in request.files.getlist("file"):
     #   print(file)
      #  filename = file.filename
       # destination = "/".join([target, filename])
        #print(destination)
        #file.save(destination)

    #return render_template("complete.html")

							
if __name__ == '__main__':
	app.run (debug=True)

