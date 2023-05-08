#!/usr/bin/env python3

"""
.. module:: convert.py
   :synopsis: simple script used to convert ipynb to html and py

"""

import subprocess, glob, os

def runNotebook ( nb : str, cmd : str ) -> None:
    """ run a single notebook """
    execute=""
    if "interactivePlots" in nb:
        execute = " --execute"
    cmd1="%s%s --to html %s" % ( cmd, execute, nb )
    cmd2="%s --to python %s" % ( cmd, nb )
    print ( "convert: %s" % cmd1 )
    subprocess.getoutput ( cmd1 )
    subprocess.getoutput ( cmd2 )

def convert( forceConvert : bool = False ) -> None:
    """ convert notebooks
    :param forceConvert: if true, then run conversion if it hasnt changed.
    """
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

    for notebook in notebooks:
        m_nb = os.stat ( notebook ).st_mtime ## last modified notebook
        htmlf = notebook.replace(".ipynb",".html")
        if not htmlf in htmls or forceConvert:
            runNotebook ( notebook, cmd )
            continue
        m_html = os.stat ( htmlf ).st_mtime ## last modified html
        if m_html < m_nb or forceConvert: ## notebook changed since?
            runNotebook ( notebook, cmd )
        else:
            print ( "%s has not changed." % notebook )

if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description="Simple script that converts the notebooks to html and python scripts")
    ap.add_argument('-f', '--forceConvert', action="store_true", 
            help='force conversion, even if notebook hasnt changed' )
    args = ap.parse_args()
    convert( forceConvert = args.forceConvert )
