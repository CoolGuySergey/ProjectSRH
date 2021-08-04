

# Description: controller for SRH Clustermapper


#------------------------------------------------------------------------

# Standard:
import sys
import os

# External:
from flask import Flask, render_template, request, session, flash, redirect, url_for, g

# Local:
# To import SRHClustermapper.py (not as another module but as smth the webtool wraps around)
# insert path to scipts at 1
# 0 is the current path (or '' in REPL), which is searched first
sys.path.insert(1, '../scripts/')
import SRHClusterMapper

#------------------------------------------------------------------------


# Configuration
MAX_CONTENT_LENGTH = 1024 * 1024
# Maximum size a request body can have is 1MB
# Requests that are larger than this are discarded with a 413 status code.
UPLOAD_PATH = 'uploads'

# Create application object
app = Flask(__name__)

# Pulls in app configuration by looking for uppercase variables:
# __name__: pass the entire module to method
# from_object() method: takes object as a parameter and pass to config
# i.e. will look for variables within the object that are defined in all caps
# currently 0 vars set up
app.config.from_object(__name__)


#------------------------------------------------------------------------


# Decorators for view function below
# Define routes for view function (Bind to URL)
# Define view with function, which does something
@app.route("/")
def index():
    error = None
    status_code = 200
    return render_template("index.html", error=error), status_code


@app.route("/results", methods=["POST"])
def results():

    UploadedFile = request.files['PathToInputAln']
    if UploadedFile.filename != '':
        UploadedFile.save("uploads/TempFileIn")
    # UploadedFile is instance of class FileStorage
    # Review and validate submission before writing to current directory
    # Secure upload handling: ignore client-provided filename before passing to  save(). 
     
    class CreateArgs:
        i = "uploads/TempFileIn"
        p = request.form['Partition']
        a = float(request.form['Alpha']

    args = CreateArgs()
    SRHClusterMapper.run(args)
    
    # TempFileIn* mop up all output pngs in the current dir (except the FileIn itself in uploads)
    maps = ["Heatmap1", "Heatmap2", "Heatmap3"]
    
    return render_template("results.html", maps=maps, Alpha=args.a, Partition=args.p, PathToInputAln=args.i)
    # Passes variable posts to main.html


#------------------------------------------------------------------------


# Start the development server using the run() method
# Switch-on error handling
if __name__ == "__main__":
    app.run(debug=True)
