

# Description: controller for SRH Clustermapper


#------------------------------------------------------------------------


# Standard:
import sys
import os
import shutil
from random import SystemRandom
from string import ascii_uppercase

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
SECRET_KEY = os.urandom(5)
# Maximum size a request body can have is 1MB
# Requests that are larger than this are discarded with a 413 status code.

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


#------------------------------------------------------------------------


@app.route("/", methods=["POST"])
def upload():
    
    UploadedFile = request.files['PathToInputAln']
    if UploadedFile.filename != '':

        session["UserID"] = ''.join(SystemRandom().choice(ascii_uppercase) for _ in range(7))
        UserDirPath = os.path.join("static/uploads", str(session["UserID"]))
        # os.path.join returns str type
        # Path: static/uploads/Òa;rO/
        
        os.mkdir(UserDirPath)
        InputPath = os.path.join(UserDirPath, "TempFileIn")
        UploadedFile.save(InputPath)
        # makes dir with path static/uploads/Òa;rO/ and saves TempFileIn
        
        class ArgsClass:
            i = InputPath
            p = request.form['Partition']
            a = float(request.form['Alpha'])

        args = ArgsClass()
        args.p = SRHClusterMapper.str2bool(args.p)
        
        SRHClusterMapper.run(args)

        return redirect(url_for("results", Partition=args.p, Alpha=args.a, UserDirPath=UserDirPath))

    else:
        error = "Please provide an alignment."
        status_code = 401
        
    return render_template("index.html", error=error), status_code


#------------------------------------------------------------------------


@app.route("/results")
def results():
    
    UserDirPath = request.args.get("UserDirPath")
    
    AllItems = os.listdir(UserDirPath)
    ValidItems = []
    for Item in AllItems:
        if Item.endswith(".png"):
            ValidItems.append(Item)
        
    ImageNames = [UserDirPath + "/" + x for x in ValidItems]

    return render_template("results.html", ImageNames=ImageNames, Partition=request.args.get("Partition"), Alpha=request.args.get("Alpha"))


#------------------------------------------------------------------------


@app.route("/clear")
def clear():
    
    # Reset session
    UserDirPath = os.path.join("static/uploads", session["UserID"])
    shutil.rmtree(UserDirPath)
    session.pop("UserID", None)
    
    return redirect(url_for("index"))


#------------------------------------------------------------------------


# Start the development server using the run() method
# Switch-on error handling
if __name__ == "__main__":
    app.run(debug=True)
