#!/usr/bin/env python3

import subprocess, glob, os

cmd=subprocess.getoutput("which jupyter-nbconvert")
if cmd=="":
    cmd=subprocess.getoutput("which jupyter" )
    if cmd == "":
        cmd=subprocess.getoutput("which ipython" ) 
    if cmd != "":
        cmd =cmd + " nbconvert"

print ( "command for conversion: %s" % cmd )

notebooks=glob.glob("*.ipynb")
htmls=glob.glob("*.html")

def run ( nb ):
    cmd1="%s --to html %s" % ( cmd, nb )
    cmd2="%s --to python %s" % ( cmd, nb )
    print ( "convert: %s" % cmd1 )
    subprocess.getoutput ( cmd1 )
    subprocess.getoutput ( cmd2 )

for notebook in notebooks:
    m_nb = os.stat ( notebook ).st_mtime ## last modified notebook
    htmlf = notebook.replace(".ipynb",".html")
    if not htmlf in htmls:
        run ( notebook )
        continue
    m_html = os.stat ( htmlf ).st_mtime ## last modified html
    if m_html < m_nb: ## notebook changed since?
        run ( notebook )
    else:
        print ( "%s has not changed." % notebook )

#
#for i in `ls *.ipynb`; do
#	$CMD --to html $i || { echo "\nERROR: cannot execute $CMD nbconvert. Maybe install jupyter-nbconvert?\n\n"; exit; }
#  $CMD --to python $i;
#done
#
