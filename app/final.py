######## Melonomics Software Development Project############
# for: QMUL MSc Bioinformatics
#Name: Melon Group (Modupeh Betts, Nadim Rahman, Maddy Rhodes, Andrew Knowles)
#Date: Submitted for assement 17/2/17



#imports python os package to enable access to the sub-modules such as os.path specifying file saving path
import os
#imports modules of flask that are required
from flask import Flask, render_template, request, url_for, redirect, send_file, send_from_directory

#this initializes the flask app

app = Flask(__name__)

APP__ROOT = os.path.dirname(os.path.abspath(__file__))

############################################################################

#this returns the home page which is templated 
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
	return redirect (url_for("R_Downloads"))



#R analysis: This imports the runner.py file which invokes routput.R (the r analysis )of our files
#
@app.route("/R_analysis/")
def R_analysis (name=None):
    import runner
    print ("R is complete")
    return redirect (url_for("R_Downloads"))




# This shows a template which consists of a button to allow users to view R analysis.
@app.route("/R_Downloads/")
def R_Downloads():
	return render_template("Rdownloads.html")


#return R pdf graphs
# This displays a pdf of R analysis on the screen.
@app.route("/return_file/")
def returnFile():
    try:
	   return send_file("/Users/modoupehbetts/Documents/Software_development/mastermelon/app/HCA.pdf", attachment_filename="HCA.pdf")
    except Exception as e:
        return str(e)



# This displays a pdf of the heatmap on screen.
@app.route("/heatmap/")
def heatmap():
	   return send_file("/Users/modoupehbetts/Documents/Software_development/mastermelon/app/heatmap.pdf", attachment_filename="heatmap.pdf")
    

#################################Routes for HTML interactive plots############################################

@app.route('/scores1/')
def scores1():
    try:
        return send_from_directory ("/Users/modoupehbetts/Documents/Software_development/mastermelon/app/", "scores1.html")
    except Exception as e:
        return redirect (url_for("error"))
@app.route('/scores2/')


def scores2():
    try:
        return send_from_directory ("/Users/modoupehbetts/Documents/Software_development/mastermelon/app/", "scores2.html")
    except Exception as e:
        return str(e)

@app.route('/scores3/')
def scores3():
    try:
        return send_from_directory ("/Users/modoupehbetts/Documents/Software_development/mastermelon/app/", "scores3.html")
    except Exception as e:
        return redirect (url_for("error"))


@app.route('/scores4/')
def scores4():
    try:
        return send_from_directory ("/Users/modoupehbetts/Documents/Software_development/mastermelon/app/", "scores4.html")
    except Exception as e:
        return str(e)


@app.route('/variance/')
def variance():
    try:    
        return send_from_directory ("/Users/modoupehbetts/Documents/Software_development/mastermelon/app/", "variance.html")
    except Exception as e:
        return str(e)

@app.route('/volcanoplot/')
def volcanoplot():
    try:
        return send_from_directory ("/Users/modoupehbetts/Documents/Software_development/mastermelon/app/", "volcanoplot.html")            
    except Exception as e:
        return str(e)

@app.route('/toptable/')
def toptable():
    return send_from_directory ("/Users/modoupehbetts/Documents/Software_development/mastermelon/app/", "toptable.html")

####################################Routes for JavaScript of HTML R Outputs##########################################

@app.route('/scores1/<path:path>')
def scores_1(path):
    try:
        return send_from_directory('/Users/modoupehbetts/Documents/Software_development/mastermelon/app/', path)
    except Exception as e:
        return str(e)

@app.route('/scores2/<path:path>')
def scores_2(path):
    try:
        return send_from_directory('/Users/modoupehbetts/Documents/Software_development/mastermelon/app/', path)
    except Exception as e:
        return str(e)

@app.route('/scores3/<path:path>')
def scores_3(path):
    try:    
        return send_from_directory('/Users/modoupehbetts/Documents/Software_development/mastermelon/app/', path)
    except Exception as e:
        return str(e)

@app.route('/scores4/<path:path>')
def scores_4(path):
    try:
        return send_from_directory('/Users/modoupehbetts/Documents/Software_development/mastermelon/app/', path)
    except Exception as e:
        return str(e)

@app.route('/toptable/<path:path>')
def top_table(path):
    try:
        return send_from_directory('/Users/modoupehbetts/Documents/Software_development/mastermelon/app/', path)
    except Exception as e:
        return str(e)

@app.route('/variance/<path:path>')
def va_riance(path):
    try:
        return send_from_directory('/Users/modoupehbetts/Documents/Software_development/mastermelon/app/', path)
    except Exception as e:
        return str(e)



@app.route('/volcanoplot/<path:path>')
def volcano_plot(path):
    try:
        return send_from_directory('/Users/modoupehbetts/Documents/Software_development/mastermelon/app/', path)
    except Exception as e:
        return str(e)



###############################################################################################
#this route ask for input files control files and save in CTRL folder
@app.route("/analyse/")
def index():
    return render_template("upload.html")



@app.route("/upload", methods=['POST'])
def upload_1():
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





@app.route("/upload_2", methods=['POST'])
def upload_2():
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

@app.route("/upload_3", methods=['POST'])
def upload_3():
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

    return redirect(url_for("Blast"))


###########################################################################

#This is an alternative method of saving files that we may try later
#need to 
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

